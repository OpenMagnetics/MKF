#include "Settings.h"
#include "Painter.h"
#include "CoilWrapper.h"
#include "json.hpp"
#include "Utils.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
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


SUITE(CoilWeb) {
    bool plot = false;
    auto settings = OpenMagnetics::Settings::GetInstance();
    TEST(Test_Coil_Json_0) {
        std::string coilString = R"({"bobbin":"Dummy","functionalDescription":[{"isolationSide":"Primary","name":"Primary","numberParallels":1,"numberTurns":23,"wire":"Dummy"}]})";

        auto coilJson = json::parse(coilString);
        auto coilWrapper(coilJson);
    }

    TEST(Test_Coil_Json_1) {
        std::string coilString = R"({"_interleavingLevel":3,"_windingOrientation":"contiguous","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":1,"numberTurns":9,"wire":"Round 0.475 - Grade 1"}]})";

        auto coilJson = json::parse(coilString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::CoilFunctionalDescription>(coilJson["functionalDescription"]);
        OpenMagnetics::CoilWrapper coil;
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
        coil.wind();

        auto section = coil.get_sections_description().value()[0];
        CHECK(!std::isnan(section.get_dimensions()[0]));
        CHECK(!std::isnan(section.get_dimensions()[1]));
    }

    TEST(Test_Coil_Json_2) {
        std::string coilString = R"({"_interleavingLevel":7,"_windingOrientation":"overlapping","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":27,"numberTurns":36,"wire":"Round 0.475 - Grade 1"}]})";
        settings->set_coil_wind_even_if_not_fit(false);

        auto coilJson = json::parse(coilString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::CoilFunctionalDescription>(coilJson["functionalDescription"]);
        OpenMagnetics::CoilWrapper coil;
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
        coil.wind();

        auto section = coil.get_sections_description().value()[0];
        CHECK(!std::isnan(section.get_dimensions()[0]));
        CHECK(!std::isnan(section.get_dimensions()[1]));
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Json_2.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST(Test_Coil_Json_3) {
        std::string coilString = R"({"_interleavingLevel":7,"_windingOrientation":"contiguous","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":88,"numberTurns":1,"wire":"Round 0.475 - Grade 1"}]})";
        settings->set_coil_delimit_and_compact(false);

        auto coilJson = json::parse(coilString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::CoilFunctionalDescription>(coilJson["functionalDescription"]);
        OpenMagnetics::CoilWrapper coil;
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
        coil.wind();

        auto section = coil.get_sections_description().value()[0];
        CHECK(!std::isnan(section.get_dimensions()[0]));
        CHECK(!std::isnan(section.get_dimensions()[1]));
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {88};
        uint8_t interleavingLevel = 7;
        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Json_3.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST(Test_Coil_Json_4) {
        std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0032500000000000003,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0002835287369864788,"coordinates":[0.0095,0,0],"height":null,"radialHeight":0.0095,"sectionsAlignment":"outer or bottom","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":2700,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.116868676970209e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000724},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 21.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000757},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"21 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.620411001469214e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000767},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 20.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000831},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"20.5 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription": null, "turnsDescription":null,"_turnsAlignment":{"Primary section 0":"spread","Secondary section 0":"spread"},"_layersOrientation":{"Primary section 0":"overlapping","Secondary section 0":"overlapping"}})";

        std::vector<size_t> pattern = {0, 1};
        std::vector<double> proportionPerWinding = {0.5, 0.5};
        size_t repetitions = 2;

        auto coilJson = json::parse(coilString);

        auto coilFunctionalDescription = std::vector<OpenMagnetics::CoilFunctionalDescription>(coilJson["functionalDescription"]);
        OpenMagnetics::CoilWrapper coil;
        coil.set_bobbin(coilJson["bobbin"]);
        coil.set_functional_description(coilFunctionalDescription);

        if (coilJson["_layersOrientation"].is_object()) {
            auto layersOrientationPerSection = std::map<std::string, OpenMagnetics::WindingOrientation>(coilJson["_layersOrientation"]);
            for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
                coil.set_layers_orientation(layerOrientation, sectionName);
            }
        }
        else if (coilJson["_layersOrientation"].is_array()) {
            coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
            if (coil.get_sections_description()) {
                auto sections = coil.get_sections_description_conduction();
                auto layersOrientationPerSection = std::vector<OpenMagnetics::WindingOrientation>(coilJson["_layersOrientation"]);
                for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                    if (sectionIndex < layersOrientationPerSection.size()) {
                        coil.set_layers_orientation(layersOrientationPerSection[sectionIndex], sections[sectionIndex].get_name());
                    }
                }
            }
        }
        else {
            OpenMagnetics::WindingOrientation layerOrientation(coilJson["_layersOrientation"]);
            coil.set_layers_orientation(layerOrientation);

        }
        if (coilJson["_turnsAlignment"].is_object()) {
            auto turnsAlignmentPerSection = std::map<std::string, OpenMagnetics::CoilAlignment>(coilJson["_turnsAlignment"]);
            for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
                coil.set_turns_alignment(turnsAlignment, sectionName);
            }
        }
        else if (coilJson["_turnsAlignment"].is_array()) {
            coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
            if (coil.get_sections_description()) {
                auto sections = coil.get_sections_description_conduction();
                auto turnsAlignmentPerSection = std::vector<OpenMagnetics::CoilAlignment>(coilJson["_turnsAlignment"]);
                for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                    if (sectionIndex < turnsAlignmentPerSection.size()) {
                        coil.set_turns_alignment(turnsAlignmentPerSection[sectionIndex], sections[sectionIndex].get_name());
                    }
                }
            }
        }
        else {
            OpenMagnetics::CoilAlignment turnsAlignment(coilJson["_turnsAlignment"]);
            coil.set_turns_alignment(turnsAlignment);
        }
 
        coil.set_bobbin(coilJson["bobbin"]);
        coil.set_functional_description(coilFunctionalDescription);
        coil.wind();
        CHECK(bool(coil.get_sections_description()));
        CHECK(bool(coil.get_layers_description()));
        CHECK(bool(coil.get_turns_description()));

    }
}

SUITE(CoilSectionsDescriptionMargins) {
    bool plot = false;
    auto settings = OpenMagnetics::Settings::GetInstance();
    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered) {
        settings->reset();
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth, 0.001);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.001;
        
        settings->set_coil_fill_sections_with_margin_tape(false);
        // settings->set_coil_wind_even_if_not_fit(false);
        // settings->set_coil_try_rewind(true);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Top) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1], 0.0001);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth, 0.0001);
        CHECK_CLOSE(marginAfterMarginFill[0], marginAfterMarginNoFill[0], 0.0001);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Top_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginBeforeMargin_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Top_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Top_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Top_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_2[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(0, marginBeforeMargin_2[1]);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_2[1], sectionDimensionsAfterMarginNoFill_2[1], 0.0001);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_0[0], marginAfterMarginNoFill_0[0], 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_2[0], marginAfterMarginNoFill_2[0], 0.0001);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_2[1] > marginAfterMarginNoFill_2[1]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_CLOSE(sectionDimensionsBeforeMargin_1[1], sectionDimensionsAfterMarginNoFill_1[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsBeforeMargin_2[1], sectionDimensionsAfterMarginNoFill_2[1], 0.0001);
        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1], 0.0001);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth, 0.0001);
        CHECK_CLOSE(marginAfterMarginFill[1], marginAfterMarginNoFill[1], 0.0001);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginBeforeMargin_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Bottom_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Bottom_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Bottom_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_2[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(0, marginBeforeMargin_2[1]);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_2[1], sectionDimensionsAfterMarginNoFill_2[1], 0.0001);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1], 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1], 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_2[1], marginAfterMarginNoFill_2[1], 0.0001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_2[0] > marginAfterMarginNoFill_2[0]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_CLOSE(sectionDimensionsBeforeMargin_1[1], sectionDimensionsAfterMarginNoFill_1[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsBeforeMargin_2[1], sectionDimensionsAfterMarginNoFill_2[1], 0.0001);
        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1], 0.0001);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth, 0.0001);
        CHECK_CLOSE(marginAfterMarginFill[1], marginAfterMarginNoFill[1], 0.0001);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginBeforeMargin_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_No_Margin_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Horizontal_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_2[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(0, marginBeforeMargin_2[1]);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1], 0.0001);
        CHECK_CLOSE(sectionDimensionsAfterMarginFill_2[1], sectionDimensionsAfterMarginNoFill_2[1], 0.0001);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.0001);
        CHECK_CLOSE(marginAfterMarginFill_0[1], marginAfterMarginNoFill_0[1], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK_CLOSE(marginAfterMarginFill_2[1], marginAfterMarginNoFill_2[1], 0.0001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_2[0] > marginAfterMarginNoFill_2[0]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_CLOSE(sectionDimensionsBeforeMargin_1[1], sectionDimensionsAfterMarginNoFill_1[1], 0.0001);
        CHECK(sectionDimensionsBeforeMargin_2[1] > sectionDimensionsAfterMarginNoFill_2[1]);
        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Inner_No_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth, 0.001);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.001;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Inner_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Inner_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Inner_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowEndingWidth = windingWindowCoordinates[0] + windingWindowDimensions[0] / 2;
        auto sectionEndingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] + coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Outer_No_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1]);
        CHECK_CLOSE(windingWindowEndingWidth, sectionEndingWidth, 0.001);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.001;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Outer_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Outer_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Outer_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[1], sectionDimensionsAfterMarginNoFill[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth, 0.001);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[1] > sectionDimensionsAfterMarginNoFill[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 25, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.001;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Horizontal_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Horizontal_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[1] > sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_1[1] > sectionDimensionsAfterMarginNoFill_1[1]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginNoFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto marginAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK_CLOSE(marginAfterMarginFill_2[1], marginAfterMarginNoFill_2[1], 0.0001);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Top_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Top_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Top_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Bottom_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Bottom_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Bottom_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1], 0.0001);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Centered_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 3, margin * 0.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Inner.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Inner_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Inner_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Outer.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Outer_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Outer_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1], 0.0001);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Top_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Inner.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Inner_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Inner_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Outer.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Outer_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Outer_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1], 0.0001);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Bottom_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Centered_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Centered_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Inner.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Inner_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Inner_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Inner_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Outer.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Outer_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Outer_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Outer_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[1], marginAfterMarginNoFill_1[1], 0.0001);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_1[0]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread) {
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         "PQ 28/20",
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);
        auto sectionDimensionsBeforeMargin = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginBeforeMargin = coil.get_sections_description_conduction()[0].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = coil.get_sections_description_conduction()[0].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin[0]);
        CHECK_EQUAL(0, marginBeforeMargin[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill[0], sectionDimensionsAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[0] > marginAfterMarginNoFill[0]);
        CHECK(marginAfterMarginFill[1] > marginAfterMarginNoFill[1]);
        CHECK(sectionDimensionsBeforeMargin[0] > sectionDimensionsAfterMarginNoFill[0]);


        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins) {
        std::vector<int64_t> numberTurns = {34, 12, 10};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.0005;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;;
        
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
        auto marginBeforeMargin_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginBeforeMargin_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        auto core = OpenMagneticsTesting::get_quick_core("PQ 28/20", json::parse("[]"), 1, "Dummy");

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Spread_Three_Different_Margins_No_Margin.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginNoFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginNoFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();
        auto bobbin = coil.resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        auto windingWindowStartingWidth = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
        auto sectionStartingWidth_0 = coil.get_sections_description_conduction()[0].get_coordinates()[0] - coil.get_sections_description_conduction()[0].get_dimensions()[0] / 2;

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        OpenMagneticsTesting::check_turns_description(coil);

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_margin().value();
        auto marginAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_margin().value();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Add_Margin_Spread_No_Filling_Then_Filling_Vertical_Spread_Three_Different_Margins.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_EQUAL(0, marginBeforeMargin_0[0]);
        CHECK_EQUAL(0, marginBeforeMargin_0[1]);
        CHECK_EQUAL(0, marginBeforeMargin_1[0]);
        CHECK_EQUAL(0, marginBeforeMargin_1[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_0[1], sectionDimensionsAfterMarginNoFill_0[1]);
        CHECK_EQUAL(sectionDimensionsAfterMarginFill_1[1], sectionDimensionsAfterMarginNoFill_1[1]);
        CHECK_CLOSE(windingWindowStartingWidth, sectionStartingWidth_0, 0.001);
        CHECK(marginAfterMarginFill_0[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_0[1] > marginAfterMarginNoFill_0[1]);
        CHECK_CLOSE(marginAfterMarginFill_1[0], marginAfterMarginNoFill_1[0], 0.0001);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_1[1]);
        CHECK(marginAfterMarginFill_1[0] > marginAfterMarginNoFill_0[0]);
        CHECK(marginAfterMarginFill_1[1] > marginAfterMarginNoFill_0[1]);
        CHECK(sectionDimensionsBeforeMargin_0[0] > sectionDimensionsAfterMarginNoFill_0[0]);
        CHECK(sectionDimensionsBeforeMargin_1[0] > sectionDimensionsAfterMarginNoFill_1[0]);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }
}

SUITE(CoilSectionsDescriptionRectangular) {
    bool plot = false;
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_Wind_By_Section_Wind_By_Consecutive_Parallels) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Test_Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced) {
        std::vector<int64_t> numberTurns = {41};
        std::vector<int64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Test_Wind_By_Section_Wind_By_Full_Turns) {
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Test_Wind_By_Section_Wind_By_Full_Parallels) {
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Test_Wind_By_Section_Wind_By_Full_Parallels_Multiwinding) {
        std::vector<int64_t> numberTurns = {2, 5};
        std::vector<int64_t> numberParallels = {7, 7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Test_Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced_Vertical) {
        std::vector<int64_t> numberTurns = {41};
        std::vector<int64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_0) {
        std::vector<int64_t> numberTurns = {9};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_1) {
        std::vector<int64_t> numberTurns = {6};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_2) {
        std::vector<int64_t> numberTurns = {5};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_3) {
        std::vector<int64_t> numberTurns = {5};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_4) {
        std::vector<int64_t> numberTurns = {91};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_5) {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_6) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {43};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 5;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random) {
        settings->set_coil_try_rewind(false);
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<int64_t> numberTurns = {std::rand() % 100 + 1L};
            std::vector<int64_t> numberParallels = {std::rand() % 100 + 1L};
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            int64_t numberPhysicalTurns = numberTurns[0] * numberParallels[0];
            uint8_t interleavingLevel = uint8_t(std::rand() % 10 + 1);
            interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::CONTIGUOUS : OpenMagnetics::WindingOrientation::OVERLAPPING;

            auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
        settings->reset();
    }

    TEST(Test_Wind_By_Section_Random_Multiwinding) {
        settings->set_coil_try_rewind(false);
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<int64_t> numberTurns;
            std::vector<int64_t> numberParallels;
            int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
            for (size_t windingIndex = 0; windingIndex < std::rand() % 10 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1L);
                numberParallels.push_back(std::rand() % 100 + 1L);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
            }
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            int64_t interleavingLevel = std::rand() % 10 + 1;
            interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::CONTIGUOUS : OpenMagnetics::WindingOrientation::OVERLAPPING;
            if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
                bobbinWidth *= numberTurns.size();
            }
            else {
                bobbinHeight *= numberTurns.size();
            }

            auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
        settings->reset();
    }

    TEST(Test_Wind_By_Section_With_Insulation_Sections) {
        std::vector<int64_t> numberTurns = {23, 42};
        std::vector<int64_t> numberParallels = {2, 1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;

        auto wires = std::vector<OpenMagnetics::WireWrapper>({OpenMagnetics::find_wire_by_name("Round 0.014 - Grade 1")});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        double voltagePeakToPeak = 400;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        auto log = coil.read_log();

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, sectionOrientation);
    }

    TEST(Test_Wind_By_Section_Pattern) {
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
}

SUITE(CoilLayersDescription) {
    bool plot = false;
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_Wind_By_Layer_Wind_One_Section_One_Layer) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_One_Section_Two_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_One_Section_One_Layer_Two_Parallels) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_One_Section_Two_Layers_Two_Parallels) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_Two_Sections_Two_Layers_Two_Parallels) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_Two_Sections_One_Layer_One_Parallel) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_Two_Sections_One_Layer_Two_Parallels) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_Two_Sections_Two_Layers_One_Parallel) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Wind_Vertical_Winding_Horizontal_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layer_Wind_Vertical_Winding_Vertical_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layer_Wind_Horizontal_Winding_Horizontal_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layer_Wind_Horizontal_Winding_Vertical_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layer_Wind_Horizontal_Winding) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil);
    }

    TEST(Test_Wind_By_Layer_Random_0) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_Wind_By_Layer_Random) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<int64_t> numberTurns = {std::rand() % 10 + 1L};
            std::vector<int64_t> numberParallels = {std::rand() % 3 + 1L};
            double wireDiameter = 0.000509;
            int64_t numberMaximumTurnsPerLayer = std::rand() % 4 + 1L;
            int64_t numberMaximumLayers = std::rand() % 3 + 1L;
            uint8_t interleavingLevel = std::rand() % 10 + 1;
            int64_t numberPhysicalTurns = numberTurns[0] * numberParallels[0];
            interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
            double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
            double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

            auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
            OpenMagneticsTesting::check_layers_description(coil);
        }
    }

    TEST(Test_Wind_By_Layer_With_Insulation_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {23, 42};
        std::vector<int64_t> numberParallels = {2, 1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;

        auto wires = std::vector<OpenMagnetics::WireWrapper>({OpenMagnetics::find_wire_by_name("Round 0.014 - Grade 1")});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        double voltagePeakToPeak = 400;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        auto log = coil.read_log();

        OpenMagneticsTesting::check_layers_description(coil);
    }
}

SUITE(CoilTurnsDescription) {
    bool plot = true;
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding) {
        settings->set_coil_wind_even_if_not_fit(false);
        srand (time(NULL));
        auto numberReallyTestedWound = std::vector<int>(2, 0);
        for (size_t testIndex = 0; testIndex < 2; ++testIndex) {
            if (testIndex == 0) {
                settings->set_coil_try_rewind(false);
            }
            else {
                settings->set_coil_try_rewind(true);
            }

            for (size_t i = 0; i < 100; ++i)
            {
                std::vector<int64_t> numberTurns;
                std::vector<int64_t> numberParallels;
                int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
                for (size_t windingIndex = 0; windingIndex < std::rand() % 2 + 1UL; ++windingIndex)
                // for (size_t windingIndex = 0; windingIndex < std::rand() % 10 + 1UL; ++windingIndex)
                {
                    int64_t numberPhysicalTurnsThisWinding = std::rand() % 300 + 1UL;
                    int64_t numberTurnsThisWinding = std::rand() % 100 + 1L;
                    int64_t numberParallelsThisWinding = std::max(1.0, std::ceil(double(numberPhysicalTurnsThisWinding) / numberTurnsThisWinding));
                    numberTurns.push_back(numberTurnsThisWinding);
                    numberParallels.push_back(numberParallelsThisWinding);
                    numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
                }
                double bobbinHeight = 0.01;
                double bobbinWidth = 0.01;
                std::vector<double> bobbinCenterCoodinates = {0.05, 0, 0};
                uint8_t interleavingLevel = std::rand() % 10 + 1;
                interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
                int windingOrientationIndex = std::rand() % 2;
                OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(windingOrientationIndex).value();

                // auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::CONTIGUOUS : OpenMagnetics::WindingOrientation::OVERLAPPING;
                if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
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

        CHECK(numberReallyTestedWound[1] > numberReallyTestedWound[0]);

        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_0) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_1) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_layers_description(coil);
        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_2) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Wind_By_Turn_Random_Multiwinding_2.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.paint_coil_sections(magnetic);
            // painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_3) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_layers_description(coil);
        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_4) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Wind_By_Turn_Random_Multiwinding_4.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_5) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Wind_By_Turn_Random_Multiwinding_4.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_6) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Random_Multiwinding_7) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Wind_By_Turn_Random_Multiwinding_4.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.paint_coil_sections(magnetic);
            // painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Rectangular_No_Bobbin) {
        settings->set_coil_wind_even_if_not_fit(false);
        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0038);
        wire.set_nominal_value_conducting_height(0.00076);
        wire.set_nominal_value_outer_width(0.004);
        wire.set_nominal_value_outer_height(0.0008);
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Different_Alignment) {
        settings->set_coil_wind_even_if_not_fit(false);
        std::vector<int64_t> numberTurns = {10, 10};
        std::vector<int64_t> numberParallels = {1, 1};
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
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
            bobbinWidth *= numberTurns.size();
            bobbinCenterCoodinates[0] += bobbinWidth / 2;
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        coil.clear();
        coil.wind();
        coil.set_turns_alignment(OpenMagnetics::CoilAlignment::INNER_OR_TOP, "winding 0 section 0");
        coil.set_turns_alignment(OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM, "winding 1 section 0");
        coil.wind_by_turns();
        coil.delimit_and_compact();
        

        OpenMagneticsTesting::check_turns_description(coil);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Different_Alignment.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_coil(coil);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            // painter.paint_coil_sections(magnetic);
            // painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        settings->reset();
    }
}

SUITE(CoilTurnsDescriptionToroidalNoCompact) { 
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = false;

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Large_Layer_Toroidal) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Large_Layer_Toroidal.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_EQUAL(1, coil.get_layers_description().value().size());
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Full_Layer_Toroidal) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {58};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Full_Layer_Toroidal.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_EQUAL(1, coil.get_layers_description().value().size());
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_Two_Layers_Toroidal) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {59};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_Two_Layers_Toroidal.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {3};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset(); 
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(180, coil.get_turns_description().value()[1].get_coordinates()[1], 0.001);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Top) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {3};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Bottom) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {3};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(357, coil.get_turns_description().value()[2].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Spread) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {3};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(60, coil.get_turns_description().value()[0].get_coordinates()[1], 0.5);
        CHECK_CLOSE(180, coil.get_turns_description().value()[1].get_coordinates()[1], 0.5);
        CHECK_CLOSE(300, coil.get_turns_description().value()[2].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Contiguous_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {3, 3};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        settings->set_coil_try_rewind(false);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Contiguous_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(90, coil.get_turns_description().value()[1].get_coordinates()[1], 0.5);
        CHECK_CLOSE(270, coil.get_turns_description().value()[4].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Overlapping_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {55, 55};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Overlapping_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_Four_Sections_One_Layer_Toroidal_Overlapping_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        // settings->set_coil_delimit_and_compact(false);
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_By_Turn_Wind_Four_Sections_One_Layer_Toroidal_Overlapping_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

}

SUITE(CoilTurnsDescriptionToroidal) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = false;

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 135);
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(182, coil.get_turns_description().value()[59].get_coordinates()[1], 1);
        CHECK_CLOSE(4.25, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(327, coil.get_turns_description().value()[101].get_coordinates()[1], 1);
        CHECK_CLOSE(5.5, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        CHECK_CLOSE(299, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Bottom) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 135);
        CHECK_CLOSE(160, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(357, coil.get_turns_description().value()[59].get_coordinates()[1], 1);
        CHECK_CLOSE(32, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(356, coil.get_turns_description().value()[101].get_coordinates()[1], 1);
        CHECK_CLOSE(60, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        CHECK_CLOSE(355, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(81, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(173, coil.get_turns_description().value()[15].get_coordinates()[1], 1);
        CHECK_CLOSE(180, coil.get_turns_description().value()[16].get_coordinates()[1], 1);
        CHECK_CLOSE(272, coil.get_turns_description().value()[31].get_coordinates()[1], 1);
        CHECK_CLOSE(327, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }
 
    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        // settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(5, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(353, coil.get_turns_description().value()[59].get_coordinates()[1], 1);
        CHECK_CLOSE(354, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(3, coil.get_turns_description().value()[18].get_coordinates()[1], 1);
        CHECK_CLOSE(3, coil.get_turns_description().value()[34].get_coordinates()[1], 1);
        CHECK_CLOSE(317, coil.get_turns_description().value()[119].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Bottom) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(12, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(117, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(123, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(221, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
        // Not clearly what this combination should do, so I check nothing
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Spread) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(117, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(123, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(243, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(42, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(147, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(43, coil.get_turns_description().value()[34].get_coordinates()[1], 1);
        CHECK_CLOSE(357, coil.get_turns_description().value()[119].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Bottom) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(42, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(147, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(44, coil.get_turns_description().value()[34].get_coordinates()[1], 1);
        CHECK_CLOSE(357, coil.get_turns_description().value()[119].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Spread) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(117, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(123, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(243, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(23, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(177, coil.get_turns_description().value()[67].get_coordinates()[1], 1);
        CHECK_CLOSE(232, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        CHECK_CLOSE(329, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Bottom) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(23, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(177, coil.get_turns_description().value()[67].get_coordinates()[1], 1);
        CHECK_CLOSE(232, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        CHECK_CLOSE(336, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(23, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(177, coil.get_turns_description().value()[67].get_coordinates()[1], 1);
        CHECK_CLOSE(232, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        CHECK_CLOSE(333, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Spread) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(117, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(123, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(243, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(123, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(243, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Bottom) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Bottom.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(12, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(115, coil.get_turns_description().value()[59].get_coordinates()[1], 1);
        CHECK_CLOSE(236, coil.get_turns_description().value()[101].get_coordinates()[1], 1);
        CHECK_CLOSE(356, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Centered) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Centered.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(7, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(109, coil.get_turns_description().value()[59].get_coordinates()[1], 1);
        CHECK_CLOSE(223, coil.get_turns_description().value()[101].get_coordinates()[1], 1);
        CHECK_CLOSE(348, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(117, coil.get_turns_description().value()[17].get_coordinates()[1], 1);
        CHECK_CLOSE(123, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(243, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Different_Wires) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 20, 20};
        std::vector<int64_t> numberParallels = {1, 5, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::WireWrapper> wires;

        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 1")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.1 - Grade 2")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Litz 225x0.04 - Grade 1 - Double Served")});


        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Different_Wires.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 180);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Different_Wires) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 20, 20};
        std::vector<int64_t> numberParallels = {1, 5, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::WireWrapper> wires;

        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 1")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.1 - Grade 2")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Litz 225x0.04 - Grade 1 - Double Served")});


        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Different_Wires.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 180);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Huge_Wire) {
        std::vector<int64_t> numberTurns = {3};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        // settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::WireWrapper> wires;

        wires.push_back({OpenMagnetics::find_wire_by_name("Litz 200x0.2 - Grade 2 - Double Served")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Huge_Wire.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 3);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Rectangular_Wire) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {11, 90};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::WireWrapper> wires;

        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 1")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Rectangular_Wire.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 101);
        // Check this one manually, checking collision between two rotated rectangles is not worth it
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Rectangular_Wire) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {6, 90};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::WireWrapper> wires;

        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 1")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Rectangular_Wire.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();

        CHECK(coil.get_turns_description().value().size() == 96);
        // Check this one manually, checking collision between two rotated rectangles is not worth it
    }
}

SUITE(CoilTurnsDescriptionToroidalMargin) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = true;

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top_Margin) {
        settings->set_coil_equalize_margins(false);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 135);
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(186, coil.get_turns_description().value()[59].get_coordinates()[1], 1);
        CHECK_CLOSE(4.25, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(175, coil.get_turns_description().value()[101].get_coordinates()[1], 1);
        CHECK_CLOSE(7, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        CHECK_CLOSE(261, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Top_Margin) {
        settings->set_coil_equalize_margins(false);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(6, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(161, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(258, coil.get_turns_description().value()[102].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Top_Margin) {
        settings->set_coil_equalize_margins(false);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(33, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(188, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(332, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Centered_Top_Margin) {
        settings->set_coil_equalize_margins(false);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(20, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(174, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(318, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Margin) {
        settings->set_coil_equalize_margins(false);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(7, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(131, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(341, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Spread_Margin) {
        settings->set_coil_equalize_margins(false);
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(7, coil.get_turns_description().value()[0].get_coordinates()[1], 1);
        CHECK_CLOSE(131, coil.get_turns_description().value()[60].get_coordinates()[1], 1);
        CHECK_CLOSE(349, coil.get_turns_description().value()[134].get_coordinates()[1], 1);
        OpenMagneticsTesting::check_turns_description(coil);
    }
}

SUITE(CoilTurnsDescriptionToroidalAdditionalCoordinates) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = false;
    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Additional_Coordinates) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        auto turns = coil.get_turns_description().value();

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Additional_Coordinates.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        for (auto turn : turns) {
            CHECK(turn.get_additional_coordinates());
        }
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread_Top_Additional_Coordinates) {
        OpenMagnetics::clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {60, 42, 33};
        std::vector<int64_t> numberParallels = {1, 1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 20/10/7";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);


        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Spread_Top_Additional_Coordinates.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        auto turns = coil.get_turns_description().value();
        for (auto turn : turns) {
            CHECK(turn.get_additional_coordinates());
            if (turn.get_additional_coordinates()) {
                auto additionalCoordinates = turn.get_additional_coordinates().value();

                for (auto additionalCoordinate : additionalCoordinates){
                    CHECK(additionalCoordinate[0] < 0);
                }
            }
        }
        OpenMagneticsTesting::check_turns_description(coil);
    }
}