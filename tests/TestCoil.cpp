#include "support/Settings.h"
#include "support/Painter.h"
#include "constructive_models/Coil.h"
#include "json.hpp"
#include "support/Utils.h"
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

using namespace MAS;
using namespace OpenMagnetics;


void process_coil_configuration(OpenMagnetics::Coil& coil, json configuration, std::optional<size_t> repetitions = std::nullopt, std::optional<std::vector<double>> proportionPerWinding = std::nullopt, std::optional<std::vector<size_t>> pattern = std::nullopt) {
    if (repetitions && proportionPerWinding && pattern) {
        if (configuration["_layersOrientation"].is_object()) {
            auto layersOrientationPerSection = std::map<std::string, WindingOrientation>(configuration["_layersOrientation"]);
            for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
                coil.set_layers_orientation(layerOrientation, sectionName);
            }
        }
        else if (configuration["_layersOrientation"].is_array()) {
            coil.wind_by_sections(proportionPerWinding.value(), pattern.value(), repetitions.value());
            if (coil.get_sections_description()) {
                auto sections = coil.get_sections_description_conduction();
                auto layersOrientationPerSection = std::vector<WindingOrientation>(configuration["_layersOrientation"]);
                for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                    if (sectionIndex < layersOrientationPerSection.size()) {
                        coil.set_layers_orientation(layersOrientationPerSection[sectionIndex], sections[sectionIndex].get_name());
                    }
                }
            }
        }
        else {
            WindingOrientation layerOrientation(configuration["_layersOrientation"]);
            coil.set_layers_orientation(layerOrientation);

        }
        if (configuration["_turnsAlignment"].is_object()) {
            auto turnsAlignmentPerSection = std::map<std::string, CoilAlignment>(configuration["_turnsAlignment"]);
            for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
                coil.set_turns_alignment(turnsAlignment, sectionName);
            }
        }
        else if (configuration["_turnsAlignment"].is_array()) {
            coil.wind_by_sections(proportionPerWinding.value(), pattern.value(), repetitions.value());
            if (coil.get_sections_description()) {
                auto sections = coil.get_sections_description_conduction();
                auto turnsAlignmentPerSection = std::vector<CoilAlignment>(configuration["_turnsAlignment"]);
                for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
                    if (sectionIndex < turnsAlignmentPerSection.size()) {
                        coil.set_turns_alignment(turnsAlignmentPerSection[sectionIndex], sections[sectionIndex].get_name());
                    }
                }
            }
        }
        else {
            CoilAlignment turnsAlignment(configuration["_turnsAlignment"]);
            coil.set_turns_alignment(turnsAlignment);
        }
    }
    else {
        if (configuration.contains("_layersOrientation")) {
            coil.set_layers_orientation(configuration["_layersOrientation"]);
        }
        if (configuration.contains("_turnsAlignment")) {
            coil.set_turns_alignment(configuration["_turnsAlignment"]);
        }
    }

    if (configuration.contains("_interleavingLevel")) {
        coil.set_interleaving_level(configuration["_interleavingLevel"]);
    }
    if (configuration.contains("_windingOrientation")) {
        coil.set_winding_orientation(configuration["_windingOrientation"]);
    }
    if (configuration.contains("_sectionAlignment")) {
        coil.set_section_alignment(configuration["_sectionAlignment"]);
    }
    if (configuration.contains("_sectionAlignment")) {
        coil.set_section_alignment(configuration["_sectionAlignment"]);
    }

    if (configuration.contains("_interlayerInsulationThickness")) {
        coil.set_interlayer_insulation(configuration["_interlayerInsulationThickness"], std::nullopt, std::nullopt, false);
    }
    if (configuration.contains("_intersectionInsulationThickness")) {
        coil.set_intersection_insulation(configuration["_intersectionInsulationThickness"], 1, std::nullopt, std::nullopt, false);
    }
}


SUITE(CoilWeb) {
    bool plot = false;
    auto settings = Settings::GetInstance();
    TEST(Test_Coil_Json_0) {
        std::string coilString = R"({"bobbin":"Dummy","functionalDescription":[{"isolationSide":"Primary","name":"Primary","numberParallels":1,"numberTurns":23,"wire":"Dummy"}]})";

        auto coilJson = json::parse(coilString);
        auto Coil(coilJson);
    }

    TEST(Test_Coil_Json_1) {
        std::string coilString = R"({"_interleavingLevel":3,"_windingOrientation":"contiguous","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":1,"numberTurns":9,"wire":"Round 0.475 - Grade 1"}]})";

        auto coilJson = json::parse(coilString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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

        CHECK(coil.get_functional_description().size() > 0);
        to_json(coilJson, coil);
        CHECK(coilJson["functionalDescription"].size() > 0);

        coil.wind();

        auto section = coil.get_sections_description().value()[0];
        CHECK(!std::isnan(section.get_dimensions()[0]));
        CHECK(!std::isnan(section.get_dimensions()[1]));
    }

    TEST(Test_Coil_Json_2) {
        std::string coilString = R"({"_interleavingLevel":7,"_windingOrientation":"overlapping","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":27,"numberTurns":36,"wire":"Round 0.475 - Grade 1"}]})";
        settings->set_coil_wind_even_if_not_fit(false);

        auto coilJson = json::parse(coilString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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
        coil.wind();

        auto section = coil.get_sections_description().value()[0];
        CHECK(!std::isnan(section.get_dimensions()[0]));
        CHECK(!std::isnan(section.get_dimensions()[1]));
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    TEST(Test_Coil_Json_3) {
        std::string coilString = R"({"_interleavingLevel":7,"_windingOrientation":"contiguous","_layersOrientation":"overlapping","_turnsAlignment":"centered","_sectionAlignment":"centered","bobbin":{"processedDescription":{"columnDepth":0.005,"columnShape":"round","columnThickness":0.001,"wallThickness":0.001,"windingWindows":[{"coordinates":[0.01,0.0,0.0],"height":0.01,"width":0.01}]}},"functionalDescription":[{"isolationSide":"primary","name":"winding 0","numberParallels":88,"numberTurns":1,"wire":"Round 0.475 - Grade 1"}]})";
        settings->set_coil_delimit_and_compact(false);

        auto coilJson = json::parse(coilString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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
        coil.wind();

        auto section = coil.get_sections_description().value()[0];
        CHECK(!std::isnan(section.get_dimensions()[0]));
        CHECK(!std::isnan(section.get_dimensions()[1]));
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {88};
        uint8_t interleavingLevel = 7;
        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    TEST(Test_Coil_Json_4) {
        std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0032500000000000003,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0002835287369864788,"coordinates":[0.0095,0,0],"height":null,"radialHeight":0.0095,"sectionsAlignment":"outer or bottom","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":2700,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.116868676970209e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000724},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 21.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000757},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"21 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":27,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":4.620411001469214e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000767},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 20.5 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000831},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"20.5 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription": null, "turnsDescription":null,"_turnsAlignment":{"Primary section 0":"spread","Secondary section 0":"spread"},"_layersOrientation":{"Primary section 0":"overlapping","Secondary section 0":"overlapping"}})";

        std::vector<size_t> pattern = {0, 1};
        std::vector<double> proportionPerWinding = {0.5, 0.5};
        size_t repetitions = 2;

        auto coilJson = json::parse(coilString);

        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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
                auto layersOrientationPerSection = std::vector<WindingOrientation>(coilJson["_layersOrientation"]);
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
                auto turnsAlignmentPerSection = std::vector<CoilAlignment>(coilJson["_turnsAlignment"]);
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
        CHECK(bool(coil.get_sections_description()));
        CHECK(bool(coil.get_layers_description()));
        CHECK(bool(coil.get_turns_description()));
    }

    TEST(Test_Coil_Json_5) {
        std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.004347500000000001,"columnShape":"round","columnThickness":0.0007975000000000005,"columnWidth":0.004347500000000001,"coordinates":[0,0,0],"pins":null,"wallThickness":0.0008723921229391407,"windingWindows":[{"angle":null,"area":0.000022027638255648275,"coordinates":[0.0059425,0,0],"height":0.006905215754121718,"radialHeight":null,"sectionsAlignment":"inner or top","sectionsOrientation":"overlapping","shape":"rectangular","width":0.0031899999999999993}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":7,"wire":{"coating":{"breakdownVoltage":2359,"grade":3,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":3.1172453105244723e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00063},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 0.63 - FIW 3","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0007279999999999999,"minimum":0.000705,"nominal":null},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"0.63 mm","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":19,"wire":{"coating":{"breakdownVoltage":2700,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":3.6637960384511227e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000683},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 21.5 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000716},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"21.5 AWG","strand":null,"type":"round"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,2.168404344971009e-19],"dimensions":[0.0007164999999999999,0.006635254],"fillingFactor":0.7263350166873861,"insulationMaterial":null,"name":"Primary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.0050765,0],"dimensions":[0.000025,0.006905215754121718],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Primary and Primary section 1 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Primary and Primary section 1","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,2.168404344971009e-19],"dimensions":[0.000716,0.006634754],"fillingFactor":0.7258281534517355,"insulationMaterial":null,"name":"Secondary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[0.3684210526315789],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,-4.999999995199816e-10],"dimensions":[0.000716,0.006470344999999999],"fillingFactor":0.6221384172443447,"insulationMaterial":null,"name":"Secondary section 0 layer 1","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[0.3157894736842105],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,-4.999999995199816e-10],"dimensions":[0.000716,0.006470344999999999],"fillingFactor":0.6221384172443447,"insulationMaterial":null,"name":"Secondary section 0 layer 2","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[0.3157894736842105],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.007249500000000002,0],"dimensions":[0.000025,0.006905215754121718],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Secondary and Secondary section 3 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Secondary and Secondary section 3","turnsAlignment":"spread","type":"insulation","windingStyle":null}],"sectionsDescription":[{"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,0],"dimensions":[0.0007164999999999999,0.006635254],"fillingFactor":0.6629541903904612,"layersAlignment":null,"orientation":"overlapping","margin":[0,0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.005076500000000001,0],"dimensions":[0.000025,0.006905215754121718],"fillingFactor":1,"layersAlignment":null,"orientation":"overlapping","margin":null,"name":"Insulation between Primary and Primary section 1","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.006163000000000002,0],"dimensions":[0.002148,0.006634754],"fillingFactor":0.5926870467921613,"layersAlignment":null,"orientation":"overlapping","margin":[0,0],"name":"Secondary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.0072495000000000025,0],"dimensions":[0.000025,0.006905215754121718],"fillingFactor":1,"layersAlignment":null,"orientation":"overlapping","margin":null,"name":"Insulation between Secondary and Secondary section 3","partialWindings":[],"type":"insulation","windingStyle":null}],"turnsDescription":[{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,0.0029593770000000004],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,0.0019729180000000006],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,0.0009864590000000003],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,2.168404344971009e-19],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,-0.0009864589999999999],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,-0.0019729179999999997],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.004705750000000001,-0.002959377],"dimensions":[0.0007164999999999999,0.0007164999999999999],"layer":"Primary section 0 layer 0","length":0.029567099259260342,"name":"Primary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,0.0029593770000000004],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,0.0019729180000000006],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,0.0009864590000000003],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,2.168404344971009e-19],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,-0.0009864589999999999],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,-0.0019729179999999997],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005447000000000002,-0.002959377],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 0","length":0.03422451036820722,"name":"Secondary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,0.0028771720000000003],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 1","length":0.03872327104814781,"name":"Secondary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,0.0017263030000000004],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 1","length":0.03872327104814781,"name":"Secondary parallel 0 turn 8","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,0.0005754340000000005],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 1","length":0.03872327104814781,"name":"Secondary parallel 0 turn 9","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,-0.0005754349999999994],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 1","length":0.03872327104814781,"name":"Secondary parallel 0 turn 10","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,-0.0017263039999999993],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 1","length":0.03872327104814781,"name":"Secondary parallel 0 turn 11","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006163000000000003,-0.0028771729999999994],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 1","length":0.03872327104814781,"name":"Secondary parallel 0 turn 12","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,0.0028771720000000003],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 2","length":0.04322203172808839,"name":"Secondary parallel 0 turn 13","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,0.0017263030000000004],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 2","length":0.04322203172808839,"name":"Secondary parallel 0 turn 14","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,0.0005754340000000005],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 2","length":0.04322203172808839,"name":"Secondary parallel 0 turn 15","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,-0.0005754349999999994],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 2","length":0.04322203172808839,"name":"Secondary parallel 0 turn 16","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,-0.0017263039999999993],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 2","length":0.04322203172808839,"name":"Secondary parallel 0 turn 17","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006879000000000002,-0.0028771729999999994],"dimensions":[0.000716,0.000716],"layer":"Secondary section 0 layer 2","length":0.04322203172808839,"name":"Secondary parallel 0 turn 18","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"}],"_turnsAlignment":{"Primary section 0":"spread","Secondary section 0":"spread"},"_layersOrientation":{"Primary section 0":"overlapping","Secondary section 0":"overlapping"}})";

        std::vector<size_t> pattern = {0, 1};
        std::vector<double> proportionPerWinding = {0.25, 0.75};
        size_t repetitions = 2;

        auto coilJson = json::parse(coilString);

        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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
                auto layersOrientationPerSection = std::vector<WindingOrientation>(coilJson["_layersOrientation"]);
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
                auto turnsAlignmentPerSection = std::vector<CoilAlignment>(coilJson["_turnsAlignment"]);
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
        CHECK(bool(coil.get_sections_description()));
        CHECK(bool(coil.get_layers_description()));
        CHECK(bool(coil.get_turns_description()));
        CHECK(bool(coil.are_sections_and_layers_fitting()));
    }

    TEST(Test_Coil_Json_6) {
        json coilJson = json::parse(R"({"_sectionsAlignment":"spread","_turnsAlignment":"centered","bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.0075,"columnShape":"rectangular","columnThickness":0.0,"columnWidth":0.0026249999999999997,"coordinates":[0.0,0.0,0.0],"pins":null,"wallThickness":0.0,"windingWindows":[{"angle":360.0,"area":0.00017203361371057708,"coordinates":[0.0074,0.0,0.0],"height":null,"radialHeight":0.0074,"sectionsAlignment":"spread","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"isolationSide":"primary","name":"primary","numberParallels":1,"numberTurns":15,"wire":{"coating":{"breakdownVoltage":2700.0,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00125},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Round 1.25 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001316},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"1.25 mm","strand":null,"type":"round"}},{"isolationSide":"secondary","name":"secondary","numberParallels":1,"numberTurns":15,"wire":{"coating":{"breakdownVoltage":2700.0,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00125},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Round 1.25 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001316},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"1.25 mm","strand":null,"type":"round"}}]})");
        size_t repetitions = 1;
        json proportionPerWindingJson = json::parse(R"([])");
        json patternJson = json::parse(R"([0,1])");
        json marginPairsJson = json::parse(R"([])");

        std::vector<std::vector<double>> marginPairs;

        for (auto elem : marginPairsJson) {
            std::vector<double> vectorElem;
            for (auto value : elem) {
                vectorElem.push_back(value);
            }
            marginPairs.push_back(vectorElem);
        }

        std::vector<double> proportionPerWinding = proportionPerWindingJson;
        std::vector<size_t> pattern = patternJson;
        std::vector<OpenMagnetics::Winding> coilFunctionalDescription;
        for (auto elem : coilJson["functionalDescription"]) {
            coilFunctionalDescription.push_back(OpenMagnetics::Winding(elem));
        }
        OpenMagnetics::Coil coil;
        coil.set_bobbin(coilJson["bobbin"]);
        coil.set_functional_description(coilFunctionalDescription);
        coil.preload_margins(marginPairs);
        if (coilJson.contains("_layersOrientation")) {

            if (coilJson["_layersOrientation"].is_object()) {
                std::map<std::string, WindingOrientation> layersOrientationPerSection;
                for (auto [key, value] : coilJson["_layersOrientation"].items()) {
                    layersOrientationPerSection[key] = value;
                }

                for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
                    coil.set_layers_orientation(layerOrientation, sectionName);
                }
            }
            else if (coilJson["_layersOrientation"].is_array()) {
                coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
                if (coil.get_sections_description()) {
                    auto sections = coil.get_sections_description_conduction();

                    std::vector<WindingOrientation> layersOrientationPerSection;
                    for (auto elem : coilJson["_layersOrientation"]) {
                        layersOrientationPerSection.push_back(WindingOrientation(elem));
                    }

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
        }

        if (coilJson.contains("_turnsAlignment")) {
            if (coilJson["_turnsAlignment"].is_object()) {
                std::map<std::string, CoilAlignment> turnsAlignmentPerSection;
                for (auto [key, value] : coilJson["_turnsAlignment"].items()) {
                    turnsAlignmentPerSection[key] = value;
                }


                for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
                    coil.set_turns_alignment(turnsAlignment, sectionName);
                }
            }
            else if (coilJson["_turnsAlignment"].is_array()) {
                coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
                if (coil.get_sections_description()) {
                    auto sections = coil.get_sections_description_conduction();

                    std::vector<CoilAlignment> turnsAlignmentPerSection;
                    for (auto elem : coilJson["_turnsAlignment"]) {
                        turnsAlignmentPerSection.push_back(CoilAlignment(elem));
                    }

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
        }

        if (proportionPerWinding.size() == coilFunctionalDescription.size()) {
            if (pattern.size() > 0 && repetitions > 0) {
                coil.wind(proportionPerWinding, pattern, repetitions);
            }
            else if (repetitions > 0) {
                coil.wind(repetitions);
            }
            else {
                coil.wind();
            }
        }
        else {
            if (pattern.size() > 0 && repetitions > 0) {
                coil.wind(pattern, repetitions);
            }
            else if (repetitions > 0) {
                coil.wind(repetitions);
            }
            else {
                coil.wind();
            }
        }


        // std::cout << bool(coil.get_sections_description()) << std::endl;
        // std::cout << bool(coil.get_layers_description()) << std::endl;
        // std::cout << bool(coil.get_turns_description()) << std::endl;
        if (!coil.get_turns_description()) {
            throw std::runtime_error("Turns not created");
        }

        json result;
        to_json(result, coil);
    }

    TEST(Test_Coil_Json_7) {
        // json coilJson = json::parse(R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00356,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0022725,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0000637587014444212,"coordinates":[0.004505,0,0],"height":null,"radialHeight":0.004505,"sectionsAlignment":"inner or top","sectionsOrientation":"overlapping","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":3,"numberTurns":55,"wire":{"coating":{"breakdownVoltage":1220,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":8.042477193189871e-8},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000323,"minimum":0.00031800000000000003,"nominal":0.00032},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 28.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000356,"minimum":0.00033800000000000003,"nominal":0.000347},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"28 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread"],"_layersOrientation":["overlapping"]})");
        json coilJson = json::parse(R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00356,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0022725,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0000637587014444212,"coordinates":[0.004505,0,0],"height":null,"radialHeight":0.004505,"sectionsAlignment":"inner or top","sectionsOrientation":"overlapping","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":3,"numberTurns":55,"wire":{"coating":{"breakdownVoltage":1220,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":8.042477193189871e-8},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000323,"minimum":0.00031800000000000003,"nominal":0.00032},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 28.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000356,"minimum":0.00033800000000000003,"nominal":0.000347},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"28 AWG","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread"],"_layersOrientation":["overlapping"]})");
        size_t repetitions = 1;
        json proportionPerWindingJson = json::parse(R"([1])");
        json patternJson = json::parse(R"([0])");
        json marginPairsJson = json::parse(R"([])");

        std::vector<std::vector<double>> marginPairs;

        for (auto elem : marginPairsJson) {
            std::vector<double> vectorElem;
            for (auto value : elem) {
                vectorElem.push_back(value);
            }
            marginPairs.push_back(vectorElem);
        }

        std::vector<double> proportionPerWinding = proportionPerWindingJson;
        std::vector<size_t> pattern = patternJson;
        std::vector<OpenMagnetics::Winding> coilFunctionalDescription;
        for (auto elem : coilJson["functionalDescription"]) {
            coilFunctionalDescription.push_back(OpenMagnetics::Winding(elem));
        }
        OpenMagnetics::Coil coil;
        coil.set_bobbin(coilJson["bobbin"]);
        coil.set_functional_description(coilFunctionalDescription);
        coil.preload_margins(marginPairs);
        if (coilJson.contains("_layersOrientation")) {

            if (coilJson["_layersOrientation"].is_object()) {
                std::map<std::string, WindingOrientation> layersOrientationPerSection;
                for (auto [key, value] : coilJson["_layersOrientation"].items()) {
                    layersOrientationPerSection[key] = value;
                }

                for (auto [sectionName, layerOrientation] : layersOrientationPerSection) {
                    coil.set_layers_orientation(layerOrientation, sectionName);
                }
            }
            else if (coilJson["_layersOrientation"].is_array()) {
                coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
                if (coil.get_sections_description()) {
                    auto sections = coil.get_sections_description_conduction();

                    std::vector<WindingOrientation> layersOrientationPerSection;
                    for (auto elem : coilJson["_layersOrientation"]) {
                        layersOrientationPerSection.push_back(WindingOrientation(elem));
                    }

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
        }

        if (coilJson.contains("_turnsAlignment")) {
            if (coilJson["_turnsAlignment"].is_object()) {
                std::map<std::string, CoilAlignment> turnsAlignmentPerSection;
                for (auto [key, value] : coilJson["_turnsAlignment"].items()) {
                    turnsAlignmentPerSection[key] = value;
                }


                for (auto [sectionName, turnsAlignment] : turnsAlignmentPerSection) {
                    coil.set_turns_alignment(turnsAlignment, sectionName);
                }
            }
            else if (coilJson["_turnsAlignment"].is_array()) {
                coil.wind_by_sections(proportionPerWinding, pattern, repetitions);
                if (coil.get_sections_description()) {
                    auto sections = coil.get_sections_description_conduction();

                    std::vector<CoilAlignment> turnsAlignmentPerSection;
                    for (auto elem : coilJson["_turnsAlignment"]) {
                        turnsAlignmentPerSection.push_back(CoilAlignment(elem));
                    }

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
        }

        if (proportionPerWinding.size() == coilFunctionalDescription.size()) {
            if (pattern.size() > 0 && repetitions > 0) {
                coil.wind(proportionPerWinding, pattern, repetitions);
            }
            else if (repetitions > 0) {
                coil.wind(repetitions);
            }
            else {
                coil.wind();
            }
        }
        else {
            if (pattern.size() > 0 && repetitions > 0) {
                coil.wind(pattern, repetitions);
            }
            else if (repetitions > 0) {
                coil.wind(repetitions);
            }
            else {
                coil.wind();
            }
        }

        if (!coil.get_turns_description()) {
            throw std::runtime_error("Turns not created");
        }

        json result;
        to_json(result, coil);
    }

    TEST(Test_Coil_Json_8) {
        json coilJson = json::parse(R"({"bobbin": {"processedDescription": {"columnDepth": 0.0037755, "columnShape": "rectangular", "columnThickness": 0.0009500000000000003, "wallThickness": 0.0008999999999999998, "windingWindows": [{"area": 4.283999999999999e-05, "coordinates": [0.0055, 0.0, 0.0], "height": 0.0126, "sectionsAlignment": "inner or top", "sectionsOrientation": "overlapping", "shape": "rectangular", "width": 0.0033999999999999994}], "columnWidth": 0.0038000000000000004, "coordinates": [0.0, 0.0, 0.0]}}, "functionalDescription": [{"isolationSide": "primary", "name": "PRI", "numberParallels": 1, "numberTurns": 10, "wire": {"type": "round", "conductingDiameter": {"maximum": 0.000257, "minimum": 0.000251, "nominal": 0.000254}, "material": "copper", "outerDiameter": {"maximum": 0.000283999999999, "minimum": 0.00026900000000000003, "nominal": 0.000277}, "coating": {"breakdownVoltage": 1190.0, "grade": 1, "material": {"dielectricStrength": [{"value": 160000000.0, "temperature": 23.0, "thickness": 2.5e-05}], "name": "ETFE", "aliases": ["Tefzel ETFE"], "composition": "Ethylene Tetrafluoroethylene", "manufacturer": "Chemours", "meltingPoint": 220.0, "relativePermittivity": 2.7, "resistivity": [{"value": 1000000000000000.0, "temperature": 170.0}], "specificHeat": 1172.0, "temperatureClass": 155.0, "thermalConductivity": 0.24}, "type": "enamelled"}, "manufacturerInfo": {"name": "Elektrisola"}, "name": "Round 30.0 - Single Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "30 AWG"}}, {"isolationSide": "secondary", "name": "SEC", "numberParallels": 1, "numberTurns": 10, "wire": {"type": "round", "conductingDiameter": {"maximum": 0.00022899999999900002, "minimum": 0.000224, "nominal": 0.00022600000000000002}, "material": "copper", "outerDiameter": {"maximum": 0.000254, "minimum": 0.00023899999999900002, "nominal": 0.000245999999999}, "coating": {"breakdownVoltage": 1020.0, "grade": 1, "material": {"dielectricStrength": [{"value": 160000000.0, "temperature": 23.0, "thickness": 2.5e-05}], "name": "ETFE", "aliases": ["Tefzel ETFE"], "composition": "Ethylene Tetrafluoroethylene", "manufacturer": "Chemours", "meltingPoint": 220.0, "relativePermittivity": 2.7, "resistivity": [{"value": 1000000000000000.0, "temperature": 170.0}], "specificHeat": 1172.0, "temperatureClass": 155.0, "thermalConductivity": 0.24}, "type": "enamelled"}, "manufacturerInfo": {"name": "Elektrisola"}, "name": "Round 31.0 - Single Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "31 AWG"}}, {"isolationSide": "primary", "name": "AUX", "numberParallels": 1, "numberTurns": 10, "wire": {"type": "round", "conductingDiameter": {"maximum": 0.000257, "minimum": 0.000251, "nominal": 0.000254}, "material": "copper", "outerDiameter": {"maximum": 0.000283999999999, "minimum": 0.00026900000000000003, "nominal": 0.000277}, "coating": {"breakdownVoltage": 1190.0, "grade": 1, "material": {"dielectricStrength":[{"value": 160000000.0, "temperature": 23.0, "thickness": 2.5e-05}], "name": "ETFE", "aliases": ["Tefzel ETFE"], "composition": "Ethylene Tetrafluoroethylene", "manufacturer": "Chemours", "meltingPoint": 220.0, "relativePermittivity": 2.7, "res sistivity": [{"value": 1000000000000000.0, "temperature": 170.0}], "specificHeat": 1172.0, "temperatureClass": 155.0, "thermalConductivity": 0.24}, "type": "enamelled"}, "manufacturerInfo": {"name": "Elektrisola"}, "name": "Round 30.0 - Single B Build", "numberConductors": 1, "standard": "NEMA MW 1000 C", "standardName": "30 AWG"}}], "layersDescription": [{"coordinates": [0.0039385, 0.0], "dimensions": [0.000277, 0.010357], "name": "PRI section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"parallelsProportion": [0.5], "winding": "PRI"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.10992063492063493, "insulationMaterial": {"dielectricStrength": [{"value": 303000000.0, "temperature": 23.0, "thickness": 2.5e-05}, {"value": 240000000.0, "temperature": 23.0, "thickness": 5e-05}, {"value": 201000000.0, "temperature": 23.0, "thickness": 7.5e-05}, {"value": 154000000.0, "temperature": 23.0, "thickness": 0.000125}], "name": "Kapton HN", "aliases": [], "composition": "Polyimide", "manufacturer": "DuPont", "meltingPoint": 520.0, "relativePermittivity": 3.4, "resistivity": [{"value": 1500000000000000.0}], "specificHeat": 1090.0, "temperatureClass": 180.0, "thermalConductivity": 0.2}, "section": "PRI section 0", "turnsAlignment": "spread", "windingStyle": "windByConsecutiveParallels"}, {"coordinates": [0.0040895, 0.0], "dimensions": [2.5e-05, 0.0126], "name": "Insulation between PRI and SEC section 1 layer 0", "orientation": "overlapping", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0, "insulationMaterial": {"dielectricStrength": [{"value": 303000000.0, "temperature": 23.0, "thickness": 2.5e-05}, {"value": 240000000.0, "temperature": 23.0, "thickness": 5e-05}, {"value": 201000000.0, "temperature": 23.0, "thickness": 7.5e-05}, {"value": 154000000.0, "temperature": 23.0, "thickness": 0.000125}], "name": "Kapton HN", "aliases": [], "composition": "Polyimide", "manufacturer": "DuPont", "meltingPoint": 520.0, "relativePermittivity": 3.4, "resistivity": [{"value": 1500000000000000.0}], "specificHeat": 1090.0, "temperatureClass": 180.0, "thermalConductivity": 0.2}, "section": "Insulation between PRI and SEC section 1", "turnsAlignment": "spread"}, {"coordinates": [0.004224999999999501, 0.0], "dimensions": [0.000245999999999, 0.011585999999999], "name": "SEC section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"parallelsProportion": [1.0], "winding": "SEC"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.1952380952373016, "insulationMaterial": {"dielectricStrength": [{"value": 303000000.0, "temperature": 23.0, "thickness": 2.5e-05}, {"value": 240000000.0, "temperature": 23.0, "thickness": 5e-05}, {"value": 201000000.0, "temperature": 23.0, "thickness": 7.5e-05}, {"value": 154000000.0, "temperature": 23.0, "thickness": 0.000125}], "name": "Kapton HN", "aliases": [], "composition": "Polyimide", "manufacturer": "DuPont", "meltingPoint": 520.0, "relativePermittivity": 3.4, "resistivity": [{"value": 1500000000000000.0}], "specificHeat": 1090.0, "temperatureClass": 180.0, "thermalConductivity": 0.2}, "section": "SEC section 0", "turnsAlignment": "spread", "windingStyle": "windByConsecutiveTurns"}, {"coordinates": [0.004360499999999001, 0.0], "dimensions": [2.5e-05, 0.0126], "name": "Insulation between SEC and PRI section 3 layer 0", "orientation": "overlapping", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0, "section": "Insulation between SEC and PRI section 3", "turnsAlignment": "spread"}, {"coordinates": [0.0045114999999990016, 0.0], "dimensions": [0.000277, 0.010357], "name": "PRI section 1 layer 0", "orientation": "overlapping", "partialWindings": [{"parallelsProportion": [0.5], "winding": "PRI"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.10992063492063493, "section": "PRI section 1", "turnsAlignment": "spread", "windingStyle": "windByConsecutiveParallels"}, {"coordinates": [0.004662499999999002, 0.0], "dimensions": [2.5e-05, 0.0126], "name": "Insulation between PRI and AUX section 5 layer 0", "orientation": "overlapping", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0, "section": "Insulation between PRI and AUX section 5", "turnsAlignment": "spread"}, {"coordinates": [0.004813499999999002, 0.0], "dimensions": [0.000277, 0.011616999999999999], "name": "AUX section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"parallelsProportion": [1.0], "winding": "AUX"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.21984126984126987, "section": "AUX section 0", "turnsAlignment": "spread", "windingStyle": "windByConsecutiveTurns"}, {"coordinates": [0.004964499999999002, 0.0], "dimensions": [2.5e-05, 0.0126], "name": "Insulation between AUX and PRI section 7 layer 0", "orientation": "overlapping", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0, "section": "Insulation between AUX and PRI section 7", "turnsAlignment": "spread"}], "sectionsDescription": [{"coordinates": [0.0039385, 0.0], "dimensions": [0.000277, 0.010357], "layersOrientation": "overlapping", "name": "PRI section 0", "partialWindings": [{"parallelsProportion": [0.5], "winding": "PRI"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.05494375499265723, "margin": [0.0, 0.0], "windingStyle": "windByConsecutiveParallels"}, {"coordinates": [0.0040895, 0.0], "dimensions": [2.5e-05, 0.0126], "layersOrientation": "overlapping", "name": "Insulation between PRI and SEC section 1", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0}, {"coordinates": [0.004224999999999501, 0.0], "dimensions": [0.00024599999999900006, 0.011585999999999], "layersOrientation": "overlapping", "name": "SEC section 0", "partialWindings": [{"parallelsProportion": [1.0], "winding": "SEC"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.0433340624416858, "margin": [0.0, 0.0], "windingStyle": "windByConsecutiveTurns"}, {"coordinates": [0.0043604999999990015, 0.0], "dimensions": [2.5e-05, 0.0126], "layersOrientation": "overlapping", "name": "Insulation between SEC and PRI section 3", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0}, {"coordinates": [0.0045114999999990016, 0.0], "dimensions": [0.000277, 0.010357], "layersOrientation": "overlapping", "name": "PRI section 1", "partialWindings": [{"parallelsProportion": [0.5], "winding": "PRI"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.056211687019914226, "margin": [0.0, 0.0], "windingStyle": "windByConsecutiveParallels"}, {"coordinates": [0.004662499999999002, 0.0], "dimensions": [2.5e-05, 0.0126], "layersOrientation": "overlapping", "name": "Insulation between PRI and AUX section 5", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0}, {"coordinates": [0.004813499999999002, 0.0], "dimensions": [0.00027699999999999996, 0.011616999999999999], "layersOrientation": "overlapping", "name": "AUX section 0", "partialWindings": [{"parallelsProportion": [1.0], "winding": "AUX"}], "type": "conduction", "coordinateSystem": "cartesian", "fillingFactor": 0.054331048199001766, "margin": [0.0, 0.0], "windingStyle": "windByConsecutiveTurns"}, {"coordinates": [0.004964499999999002, 0.0], "dimensions": [2.5e-05, 0.0126], "layersOrientation": "overlapping", "name": "Insulation between AUX and PRI section 7", "partialWindings": [], "type": "insulation", "coordinateSystem": "cartesian", "fillingFactor": 1.0}], "turnsDescription": [{"coordinates": [0.0039385, 0.00504], "length": 0.031172221165044374, "name": "PRI parallel 0 turn 0", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 0"}, {"coordinates": [0.0039385, 0.00252], "length": 0.031172221165044374, "name": "PRI parallel 0 turn 1", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 0"}, {"coordinates": [0.0039385, 0.0], "length": 0.031172221165044374, "name": "PRI parallel 0 turn 2", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 0"}, {"coordinates": [0.0039385, -0.00252], "length": 0.031172221165044374, "name": "PRI parallel 0 turn 3", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 0"}, {"coordinates": [0.0039385, -0.00504], "length": 0.031172221165044374, "name": "PRI parallel 0 turn 4", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 0"}, {"coordinates": [0.004224999999999501, 0.00567], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 0", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, 0.00441], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 1", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, 0.00315], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 2", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, 0.00189], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 3", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, 0.0006299999999999999], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 4", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, -0.0006300000000000001], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 5", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, -0.0018900000000000002], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 6", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, -0.00315], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 7", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, -0.00441], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 8", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.004224999999999501, -0.00567], "length": 0.03297235375554819, "name": "SEC parallel 0 turn 9", "parallel": 0, "winding": "SEC", "coordinateSystem": "cartesian", "dimensions": [0.000245999999999, 0.000245999999999], "layer": "SEC section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "SEC section 0"}, {"coordinates": [0.0045114999999990016, 0.00504], "length": 0.034772486346052005, "name": "PRI parallel 0 turn 5", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 1 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 1"}, {"coordinates": [0.0045114999999990016, 0.00252], "length": 0.034772486346052005, "name": "PRI parallel 0 turn 6", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 1 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 1"}, {"coordinates": [0.0045114999999990016, 0.0], "length": 0.034772486346052005, "name": "PRI parallel 0 turn 7", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 1 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 1"}, {"coordinates": [0.0045114999999990016, -0.00252], "length": 0.034772486346052005, "name": "PRI parallel 0 turn 8", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 1 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 1"}, {"coordinates": [0.0045114999999990016, -0.00504], "length": 0.034772486346052005, "name": "PRI parallel 0 turn 9", "parallel": 0, "winding": "PRI", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "PRI section 1 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "PRI section 1"}, {"coordinates": [0.004813499999999002, 0.00567], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 0", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, 0.00441], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 1", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, 0.00315], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 2", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, 0.00189], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 3", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, 0.0006299999999999999], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 4", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, -0.0006300000000000001], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 5", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, -0.0018900000000000002], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 6", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, -0.00315], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 7", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, -0.00441], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 8", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}, {"coordinates": [0.004813499999999002, -0.00567], "length": 0.03667000830882024, "name": "AUX parallel 0 turn 9", "parallel": 0, "winding": "AUX", "coordinateSystem": "cartesian", "dimensions": [0.000277, 0.000277], "layer": "AUX section 0 layer 0", "orientation": "clockwise", "rotation": 0.0, "section": "AUX section 0"}]})");
        OpenMagnetics::Coil coil(coilJson, false);
        auto layers = coil.get_layers_description().value();

        for (auto layer : layers) {
            if (layer.get_type() == ElectricalType::INSULATION) {
                std::cout << layer.get_name() << std::endl;
                std::cout << "Insulation between SEC and PRI section 3 layer 0" << std::endl;
                std::cout << (layer.get_name() == "Insulation between SEC and PRI section 3 layer 0") << std::endl;

                auto material = OpenMagnetics::Coil::resolve_insulation_layer_insulation_material(coil, layer.get_name());
                json mierda;
                to_json(mierda, material);
                std::cout << mierda << std::endl;
            }

        }
    }

    TEST(Test_Coil_Json_9) {
        std::string coilString = R"({"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.01295, "columnShape": "round", "columnThickness": 0.0, "columnWidth": 0.01295, "coordinates": [0.0, 0.0, 0.0 ], "pins": null, "wallThickness": 0.0, "windingWindows": [{"angle": null, "area": 0.0001596, "coordinates": [0.0196, 0.0, 0.0 ], "height": 0.012, "radialHeight": null, "sectionsAlignment": "centered", "sectionsOrientation": "contiguous", "shape": "rectangular", "width": 0.0133 } ] } }, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 3, "numberTurns": 12, "wire": {"coating": {"breakdownVoltage": null, "grade": null, "material": {"aliases": ["Tefzel ETFE" ], "composition": "Ethylene Tetrafluoroethylene", "dielectricStrength": [{"humidity": null, "temperature": 23.0, "thickness": 2.5e-05, "value": 160000000.0 } ], "manufacturer": "Chemours", "meltingPoint": 220.0, "name": "ETFE", "relativePermittivity": 2.7, "resistivity": [{"temperature": 170.0, "value": 1000000000000000.0 } ], "specificHeat": 1172.0, "temperatureClass": 155.0, "thermalConductivity": 0.24 }, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare" }, "conductingArea": null, "conductingDiameter": null, "conductingHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00020999999999999998 }, "conductingWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002 }, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": null, "numberConductors": 1, "outerDiameter": null, "outerHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00021020999999999995 }, "outerWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002002 }, "standard": null, "standardName": null, "strand": null, "type": "rectangular" } }, {"connections": null, "isolationSide": "secondary", "name": "Secondary", "numberParallels": 3, "numberTurns": 15, "wire": {"coating": {"breakdownVoltage": null, "grade": null, "material": {"aliases": [], "composition": "Polyurethane", "dielectricStrength": [{"humidity": null, "temperature": null, "thickness": 0.0001, "value": 25000000.0 } ], "manufacturer": "MWS Wire", "meltingPoint": null, "name": "Polyurethane 155", "relativePermittivity": 3.7, "resistivity": [{"temperature": null, "value": 1e+16 } ], "specificHeat": null, "temperatureClass": 155.0, "thermalConductivity": null }, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled" }, "conductingArea": null, "conductingDiameter": null, "conductingHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00020999999999999998 }, "conductingWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002 }, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": null, "numberConductors": 1, "outerDiameter": null, "outerHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00021020999999999995 }, "outerWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002002 }, "standard": null, "standardName": null, "strand": null, "type": "rectangular" } } ], "layersOrientation": "contiguous", "turnsAlignment": "spread" })";
        auto coilJson = json::parse(coilString);
        size_t repetitions = 1;
        double insulationThickness = 0.10 / 1000;
        std::string proportionPerWindingString = "[]";
        std::string patternString = "[0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]";

        std::vector<double> proportionPerWinding = json::parse(proportionPerWindingString);
        std::vector<size_t> pattern = json::parse(patternString);
        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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

    TEST(Test_Coil_Json_10) {
        std::string coilString = R"({"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.0047, "columnShape": "round", "columnThickness": 0.0, "columnWidth": 0.0047, "coordinates": [0.0, 0.0, 0.0], "pins": null, "wallThickness": 0.0, "windingWindows": [{"angle": null, "area": 3.813000000000001e-05, "coordinates": [0.007775000000000001, 0.0, 0.0], "height": 0.006200000000000001, "radialHeight": null, "sectionsAlignment": "centered", "sectionsOrientation": "contiguous", "shape": "rectangular", "width": 0.006150000000000001}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 33, "wire": {"coating": {"breakdownVoltage": null, "grade": null, "material": {"aliases": ["Tefzel ETFE"], "composition": "Ethylene Tetrafluoroethylene", "dielectricStrength": [{"humidity": null, "temperature": 23.0, "thickness": 2.5e-05, "value": 160000000.0}], "manufacturer": "Chemours", "meltingPoint": 220.0, "name": "ETFE", "relativePermittivity": 2.7, "resistivity": [{"temperature": 170.0, "value": 1000000000000000.0}], "specificHeat": 1172.0, "temperatureClass": 155.0, "thermalConductivity": 0.24}, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare"}, "conductingArea": null, "conductingDiameter": null, "conductingHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00013900000000000002}, "conductingWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.003705}, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": null, "numberConductors": 1, "outerDiameter": null, "outerHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000139139}, "outerWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.0037087049999999996}, "standard": null, "standardName": null, "strand": null, "type": "rectangular"}}, {"connections": null, "isolationSide": "secondary", "name": "Secondary", "numberParallels": 1, "numberTurns": 30, "wire": {"coating": {"breakdownVoltage": null, "grade": null, "material": {"aliases": [], "composition": "Polyurethane", "dielectricStrength": [{"humidity": null, "temperature": null, "thickness": 0.0001, "value": 25000000.0}], "manufacturer": "MWS Wire", "meltingPoint": null, "name": "Polyurethane 155", "relativePermittivity": 3.7, "resistivity": [{"temperature": null, "value": 1e+16}], "specificHeat": null, "temperatureClass": 155.0, "thermalConductivity": null}, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": null, "conductingHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00013900000000000002}, "conductingWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.003705}, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": null, "numberConductors": 1, "outerDiameter": null, "outerHeight": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000139139}, "outerWidth": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.0037087049999999996}, "standard": null, "standardName": null, "strand": null, "type": "rectangular"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.002984091, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 0 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.002663182, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 0 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.002367273, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 1 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 1", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.002071364, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 1 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 1", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.001775455, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 2 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 2", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.001479546, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 2 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 2", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.001183637, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 3 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 3", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.000887728, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 3 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 3", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.000591819, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 4 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 4", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 0.00029591, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 4 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 4", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, 1e-09, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 5 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 5", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.000295908, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 5 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 5", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.000591817, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 6 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 6", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.000887726, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 6 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 6", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.001183635, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 7 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 7", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.001479544, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 7 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 7", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.001775453, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 8 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 8", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.002071362, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 8 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 8", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.002367271, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 9 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 9", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.00266318, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Secondary section 9 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "section": "Secondary section 9", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.007775, -0.002984089, 0.0], "dimensions": [0.006150000000000001, 0.000139139], "fillingFactor": 1.8091243902439018, "insulationMaterial": null, "name": "Primary section 10 layer 0", "orientation": "contiguous", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "section": "Primary section 10", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0029840910000000003, 0.0], "dimensions": [0.006150000000000001, 0.00023181799999999998], "fillingFactor": 1.0858507904224277, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0028181820000000002, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 1", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0026631820000000005, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0025081820000000003, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 3", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0023672730000000005, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 1", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0022263640000000006, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 5", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.002071364000000001, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 1", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.001916364000000001, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 7", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0017754550000000009, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 2", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0016345460000000008, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 9", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0014795460000000006, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 2", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0013245460000000009, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 11", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0011836370000000008, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 3", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0010427280000000008, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 13", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0008877280000000007, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 3", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0007327280000000006, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 15", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0005918190000000006, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 4", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.00045091000000000065, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 17", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0002959100000000007, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 4", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 0.0001409100000000007, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 19", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, 1.0000000007069241e-09, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 5", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0001409079999999993, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 21", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0002959079999999993, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 5", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0004509079999999993, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 23", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0005918169999999994, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 6", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0007327259999999994, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 25", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0008877259999999994, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 6", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0010427259999999992, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 27", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0011836349999999993, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 7", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0013245439999999993, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 29", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0014795439999999995, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 7", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0016345439999999993, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 31", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0017754529999999993, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 8", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0019163619999999994, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 33", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0020713619999999993, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 8", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0022263619999999995, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 35", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.0023672709999999994, 0.0], "dimensions": [0.006150000000000001, 0.00018181799999999999], "fillingFactor": 1.3844600563978613, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 9", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.002508179999999999, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Secondary section 37", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.002663179999999999, 0.0], "dimensions": [0.006150000000000001, 0.00020999999999999998], "fillingFactor": 1.1986655168292681, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Secondary section 9", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.002818179999999999, 0.0], "dimensions": [0.006150000000000001, 0.0001], "fillingFactor": 1.0, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Primary section 39", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.007775000000000001, -0.002984088999999999, 0.0], "dimensions": [0.006150000000000001, 0.00023181799999999998], "fillingFactor": 1.0858507904224277, "group": "Default group", "layersAlignment": null, "layersOrientation": "contiguous", "margin": [0.0, 0.0], "name": "Primary section 10", "partialWindings": [{"connections": null, "parallelsProportion": [0.09090909090909091], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}], "turnsDescription": null})";
        auto coilJson = json::parse(coilString);

        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
        auto coilSectionsDescription = std::vector<Section>(coilJson["sectionsDescription"]);
        auto coilLayersDescription = std::vector<Layer>(coilJson["layersDescription"]);
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
        std::cout << bool(coil.get_turns_description()) << std::endl;
    }

    TEST(Test_Coil_Json_11) {
        std::string coilString = R"({"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.006175,"columnShape":"round","columnThickness":0,"columnWidth":0.006175,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":null,"area":0.000041283000000000004,"coordinates":[0.0098875,0,0],"height":0.00556,"radialHeight":null,"sectionsAlignment":"inner or top","sectionsOrientation":"contiguous","shape":"rectangular","width":0.007425000000000001}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"primary","numberParallels":1,"numberTurns":5,"wire":{"coating":null,"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.2293100000000003e-7},"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000522},"conductingWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0023550000000000008},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 52.20 µm","numberConductors":1,"outerDiameter":null,"outerHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000522},"outerWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0023550000000000008},"standard":"IPC-6012","standardName":"1.5 oz.","strand":null,"type":"planar"}},{"connections":null,"isolationSide":"secondary","name":"SECONDARY","numberParallels":1,"numberTurns":3,"wire":{"coating":null,"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.2449700000000003e-7},"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"conductingWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.003577500000000001},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 34.80 µm","numberConductors":1,"outerDiameter":null,"outerHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"outerWidth":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.003577500000000001},"standard":"IPC-6012","standardName":"1 oz.","strand":null,"type":"planar"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null,"_turnsAlignment":["spread","spread","spread","spread"],"_layersOrientation":["contiguous","contiguous","contiguous","contiguous"],"_interlayerInsulationThickness":0,"_intersectionInsulationThickness":0.0001})";

        std::vector<size_t> pattern = {0, 1, 0, 1};
        std::vector<double> proportionPerWinding = {0.5, 0.5};
        size_t repetitions = 1;

        auto coilJson = json::parse(coilString);

        auto coilFunctionalDescription = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
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
                auto layersOrientationPerSection = std::vector<WindingOrientation>(coilJson["_layersOrientation"]);
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
                auto turnsAlignmentPerSection = std::vector<CoilAlignment>(coilJson["_turnsAlignment"]);
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
        CHECK(bool(coil.get_sections_description()));
        CHECK(bool(coil.get_layers_description()));
        CHECK(bool(coil.get_turns_description()));
    }

    TEST(Test_Coil_Json_12) {
        std::string coilString = R"({"bobbin":"Dummy","functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":{"coating":null,"conductingArea":null,"conductingDiameter":null,"conductingHeight":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0000348},"conductingWidth":{"nominal":0.002},"edgeRadius":null,"manufacturerInfo":null,"material":"copper","name":"Planar 34.80 µm","numberConductors":1,"outerDiameter":null,"outerHeight":{"nominal":0.0000348},"outerWidth":{"nominal":0.002},"standard":"IPC-6012","standardName":"1 oz.","strand":null,"type":"planar"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null})";

        settings->set_coil_wind_even_if_not_fit(true);
        auto coilJson = json::parse(coilString);
        auto coil = OpenMagnetics::Coil(coilJson, false);
        std::vector<size_t> stackUp = {0};

        coil.set_strict(false);
        coil.wind_planar(stackUp);

        if (!coil.get_turns_description()) {
            throw std::runtime_error("Turns not created");
        }

        json result;
        to_json(result, coil);
    }
}

SUITE(CoilSectionsDescriptionMargins) {
    bool plot = true;
    auto settings = Settings::GetInstance();
    TEST(Test_Add_Margin_Centered_No_Filling_Then_Filling_Horizontal_Centered) {
        settings->reset();
        std::vector<int64_t> numberTurns = {47};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        double margin = 0.002;
        
        settings->set_coil_fill_sections_with_margin_tape(false);

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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
        auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
        auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 2, 0});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto sectionDimensionsAfterMarginFill_2 = coil.get_sections_description_conduction()[2].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
        auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
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
            Painter painter(outFile);
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 3.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);
        auto marginAfterMarginFill_2 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[2]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 2.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginNoFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginNoFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        auto sectionDimensionsAfterMarginFill = coil.get_sections_description_conduction()[0].get_dimensions();
        auto marginAfterMarginFill = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_wind_even_if_not_fit(true);
        settings->set_coil_fill_sections_with_margin_tape(false);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

        settings->set_coil_fill_sections_with_margin_tape(true);
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 2.5, margin * 2.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 0.5});
        auto sectionDimensionsAfterMarginFill_0 = coil.get_sections_description_conduction()[0].get_dimensions();
        auto sectionDimensionsAfterMarginFill_1 = coil.get_sections_description_conduction()[1].get_dimensions();
        auto marginAfterMarginFill_0 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[0]);
        auto marginAfterMarginFill_1 = OpenMagnetics::Coil::resolve_margin(coil.get_sections_description_conduction()[1]);

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
    auto settings = Settings::GetInstance();

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

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_0) {
        std::vector<int64_t> numberTurns = {9};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_1) {
        std::vector<int64_t> numberTurns = {6};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_2) {
        std::vector<int64_t> numberTurns = {5};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_3) {
        std::vector<int64_t> numberTurns = {5};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_4) {
        std::vector<int64_t> numberTurns = {91};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_5) {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
    }

    TEST(Test_Wind_By_Section_Random_6) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {43};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 5;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, WindingOrientation::CONTIGUOUS);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, WindingOrientation::CONTIGUOUS);
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
            auto windingOrientation = std::rand() % 2? WindingOrientation::CONTIGUOUS : WindingOrientation::OVERLAPPING;

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
            auto windingOrientation = std::rand() % 2? WindingOrientation::CONTIGUOUS : WindingOrientation::OVERLAPPING;
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
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
    auto settings = Settings::GetInstance();

    TEST(Test_Wind_By_Layers_Wind_One_Section_One_Layer) {
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

    TEST(Test_Wind_By_Layers_Wind_One_Section_Two_Layers) {
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

    TEST(Test_Wind_By_Layers_Wind_One_Section_One_Layer_Two_Parallels) {
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

    TEST(Test_Wind_By_Layers_Wind_One_Section_Two_Layers_Two_Parallels) {
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

    TEST(Test_Wind_By_Layers_Wind_Two_Sections_Two_Layers_Two_Parallels) {
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

    TEST(Test_Wind_By_Layers_Wind_Two_Sections_One_Layer_One_Parallel) {
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

    TEST(Test_Wind_By_Layers_Wind_Two_Sections_One_Layer_Two_Parallels) {
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

    TEST(Test_Wind_By_Layers_Wind_Two_Sections_Two_Layers_One_Parallel) {
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

    TEST(Test_Wind_By_Layers_Wind_Vertical_Winding_Horizontal_Layers) {
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

        auto windingOrientation = WindingOrientation::CONTIGUOUS;
        auto layersOrientation = WindingOrientation::CONTIGUOUS;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layers_Wind_Vertical_Winding_Vertical_Layers) {
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

        auto windingOrientation = WindingOrientation::CONTIGUOUS;
        auto layersOrientation = WindingOrientation::CONTIGUOUS;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layers_Wind_Horizontal_Winding_Horizontal_Layers) {
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

        auto windingOrientation = WindingOrientation::CONTIGUOUS;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layers_Wind_Horizontal_Winding_Vertical_Layers) {
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

        auto windingOrientation = WindingOrientation::CONTIGUOUS;
        auto layersOrientation = WindingOrientation::CONTIGUOUS;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Test_Wind_By_Layers_Wind_Horizontal_Winding) {
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

        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil);
    }

    TEST(Test_Wind_By_Layers_Random_0) {
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

    TEST(Test_Wind_By_Layers_Random) {
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

    TEST(Test_Wind_By_Layers_With_Insulation_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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

    TEST(Test_External_Insulation_Layers) {

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
}

SUITE(CoilTurnsDescription) {
    bool plot = true;
    auto settings = Settings::GetInstance();

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
                WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(windingOrientationIndex).value();

                // auto windingOrientation = std::rand() % 2? WindingOrientation::CONTIGUOUS : WindingOrientation::OVERLAPPING;
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
        auto windingOrientation = WindingOrientation::CONTIGUOUS;
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Rectangular_No_Bobbin) {
        settings->set_coil_wind_even_if_not_fit(false);
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
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers) {
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Two_Times) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0.0002);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0.0001);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0.0005);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Toroidal_Core_Contiguous) {
        // settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
        settings->reset();
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Primary) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0.0001, std::nullopt, "winding 0");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_Only_Secondary) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0.0001, std::nullopt, "winding 1");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_All_Layers_Contiguous_Layers) {
        // settings->set_coil_delimit_and_compact(false);
        // settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_interlayer_insulation(0.0001);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_intersection_insulation(0.0002, 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Interleaved) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_intersection_insulation(0.0001, 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterSections_All_Sections_Contiguous) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_intersection_insulation(0.0002, 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_intersection_insulation(0.0001, 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterSections_All_Layers_Toroidal_Core_Contiguous) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

        coil.set_intersection_insulation(0.0001, 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turn_Change_Insulation_InterLayers_And_InterSections_All_Sections) {
        settings->set_coil_wind_even_if_not_fit(true);
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
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }

    std::cout << "Mierdooooooooooooooon" << std::endl;
        coil.set_interlayer_insulation(0.00005);
        coil.set_intersection_insulation(0.0002, 1);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }
}

SUITE(CoilTurnsDescriptionToroidalNoCompact) { 
    auto settings = Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = false;

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Large_Layer_Toroidal) {
        clear_databases();
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_EQUAL(1U, coil.get_layers_description().value().size());
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Full_Layer_Toroidal) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_EQUAL(1U, coil.get_layers_description().value().size());
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_Two_Layers_Toroidal) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Centered) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset(); 
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(180, coil.get_turns_description().value()[1].get_coordinates()[1], 0.001);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Top) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(3, coil.get_turns_description().value()[0].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Bottom) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(357, coil.get_turns_description().value()[2].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_One_Section_One_Layer_Toroidal_Contiguous_Spread) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
        clear_databases();
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK_CLOSE(90, coil.get_turns_description().value()[1].get_coordinates()[1], 0.5);
        CHECK_CLOSE(270, coil.get_turns_description().value()[4].get_coordinates()[1], 0.5);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_Two_Sections_One_Layer_Toroidal_Overlapping_Centered) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        settings->set_coil_delimit_and_compact(false);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_By_Turn_Wind_Four_Sections_One_Layer_Toroidal_Overlapping_Centered) {
        clear_databases();
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }
}

SUITE(CoilTurnsDescriptionToroidal) {
    auto settings = Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = true;

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top) {
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
        // Not clearly what this combination should do, so I check nothing
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Top_Spread) {
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Bottom_Spread) {
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 180);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Different_Wires) {
        clear_databases();
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
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::Wire> wires;

        wires.push_back({find_wire_by_name("Litz 200x0.2 - Grade 2 - Double Served")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        clear_databases();
        settings->set_use_toroidal_cores(true);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 3);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Rectangular_Wire) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();
        CHECK(coil.get_turns_description().value().size() == 101);
        // Check this one manually, checking collision between two rotated rectangles is not worth it
    }

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Rectangular_Wire) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
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
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
        coil.convert_turns_to_polar_coordinates();

        CHECK(coil.get_turns_description().value().size() == 96);
        // Check this one manually, checking collision between two rotated rectangles is not worth it
    }
}

SUITE(CoilTurnsDescriptionToroidalMargin) {
    auto settings = Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = true;

    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Overlapping_Top_Margin) {
        settings->set_coil_equalize_margins(false);
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
        clear_databases();
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
    auto settings = Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    bool plot = false;
    TEST(Test_Wind_Three_Sections_Two_Layer_Toroidal_Contiguous_Spread_Top_Additional_Coordinates) {
        clear_databases();
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
        clear_databases();
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

SUITE(PlanarCoil) {
    bool plot = true;
    auto settings = Settings::GetInstance();

    TEST(Test_Wind_By_Layers_Planar_One_Layer) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(layersDescription.size(), 1);
    }

    TEST(Test_Wind_By_Layers_Planar_Two_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(layersDescription.size(), 3);
    }

    TEST(Test_Wind_By_Layers_Planar_Two_Windings) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(3U, layersDescription.size());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1, layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[2].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1, layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
    }

    TEST(Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_No_Interleaved) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(7U, layersDescription.size());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(0.5, layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[2].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(0.5, layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[4].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(0.5, layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[6].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(0.5, layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
    }

    TEST(Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_No_Interleaved_Odd_Turns) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(7U, layersDescription.size());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(2.0 / 3, layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[2].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1.0 / 3, layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[4].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(2.0 / 3, layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[6].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1.0 / 3, layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
    }

    TEST(Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(7U, layersDescription.size());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(2.0 / 3, layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[2].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(2.0 / 3, layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[4].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1.0 / 3, layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings().size());
        CHECK_EQUAL("SECONDARY", layersDescription[6].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1.0 / 3, layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
    }

    TEST(Test_Wind_By_Layers_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns_With_Insulation) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK_EQUAL(7U, layersDescription.size());
        CHECK(MAS::ElectricalType::CONDUCTION == layersDescription[0].get_type());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings().size());
        CHECK_EQUAL("PRIMARY", layersDescription[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[0].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(2.0 / 3, layersDescription[0].get_partial_windings()[0].get_parallels_proportion()[0]);

        CHECK(MAS::ElectricalType::INSULATION == layersDescription[1].get_type());

        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings().size());
        CHECK(MAS::ElectricalType::CONDUCTION == layersDescription[2].get_type());
        CHECK_EQUAL("SECONDARY", layersDescription[2].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[2].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(2.0 / 3, layersDescription[2].get_partial_windings()[0].get_parallels_proportion()[0]);

        CHECK(MAS::ElectricalType::INSULATION == layersDescription[3].get_type());

        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings().size());
        CHECK(MAS::ElectricalType::CONDUCTION == layersDescription[4].get_type());
        CHECK_EQUAL("PRIMARY", layersDescription[4].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[4].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1.0 / 3, layersDescription[4].get_partial_windings()[0].get_parallels_proportion()[0]);

        CHECK(MAS::ElectricalType::INSULATION == layersDescription[5].get_type());

        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings().size());
        CHECK(MAS::ElectricalType::CONDUCTION == layersDescription[6].get_type());
        CHECK_EQUAL("SECONDARY", layersDescription[6].get_partial_windings()[0].get_winding());
        CHECK_EQUAL(1U, layersDescription[6].get_partial_windings()[0].get_parallels_proportion().size());
        CHECK_EQUAL(1.0 / 3, layersDescription[6].get_partial_windings()[0].get_parallels_proportion()[0]);
    }

    TEST(Test_Wind_By_Turns_Planar_One_Layer) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK(coil.get_turns_description());
        auto turnsDescription = coil.get_turns_description().value();
        CHECK_EQUAL(turnsDescription.size(), 7);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    TEST(Test_Wind_By_Turns_Planar_Two_Windings_Two_Layers_Interleaved_Odd_Turns_With_Insulation) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK(coil.get_turns_description());
        if (coil.get_turns_description()) {
            auto turnsDescription = coil.get_turns_description().value();
            CHECK_EQUAL(turnsDescription.size(), 25);
            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    TEST(Test_Wind_By_Turns_Planar_Many_Layers) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK(coil.get_turns_description());
        if (coil.get_turns_description()) {
            auto turnsDescription = coil.get_turns_description().value();
            CHECK_EQUAL(turnsDescription.size(), 100);
            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    TEST(Test_Wind_By_Turns_Planar_One_Layer_Distance_To_Core) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK(coil.get_turns_description());
        auto turnsDescription = coil.get_turns_description().value();
        CHECK_EQUAL(turnsDescription.size(), 7);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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

    TEST(Test_Wind_By_Turns_Planar_Many_Layers_Magnetic_Field) {
        settings->set_coil_wind_even_if_not_fit(false);
        settings->set_coil_try_rewind(false);

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
        CHECK(coil.get_turns_description());
        if (coil.get_turns_description()) {
            auto turnsDescription = coil.get_turns_description().value();
            CHECK_EQUAL(turnsDescription.size(), 100);
            if (plot) {
                double voltagePeakToPeak = 2000;
                auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, {double(numberTurns[0]) / numberTurns[1]});
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
}

SUITE(CoilTools) {
    TEST(Test_Get_Round_Wire_From_Dc_Resistance) {
        clear_databases();
        settings->set_use_toroidal_cores(true);
        std::vector<int64_t> numberTurns = {1, 60};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "EE5";
        std::string coreMaterial = "3C97"; 
        auto emptyGapping = json::array();
        // settings->set_coil_delimit_and_compact(false);
        settings->set_coil_try_rewind(false);
        settings->set_coil_wind_even_if_not_fit(true);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        std::vector<double> dcResistances = {0.00075, 1.75};
        auto wires = coil.guess_round_wire_from_dc_resistance(dcResistances, 0.01);
        CHECK_EQUAL(wires[0].get_name().value(), "Round 0.63 - Grade 1");
        CHECK_EQUAL(wires[1].get_name().value(), "Round 0.106 - Grade 1");
        for (auto wire : wires) {
            std::cout << wire.get_name().value() << std::endl;
        }
    }
}

SUITE(CoilWindingGroups) {
    TEST(Test_Wind_By_Sections_Two_Windings_Together) {
        settings->set_coil_wind_even_if_not_fit(true);
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
        CHECK_EQUAL(1, coil.get_sections_description()->size());
        CHECK_EQUAL(2, coil.get_sections_description().value()[0].get_partial_windings().size());
        CHECK_EQUAL("winding 0", coil.get_sections_description().value()[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL("winding 1", coil.get_sections_description().value()[0].get_partial_windings()[1].get_winding());
        auto virtualFunctionalDescription = coil.virtualize_functional_description();
        CHECK_EQUAL(1, virtualFunctionalDescription.size());
        CHECK_EQUAL(numberTurns[0] + numberTurns[1], virtualFunctionalDescription[0].get_number_turns());
        CHECK_EQUAL(numberParallels[0], virtualFunctionalDescription[0].get_number_parallels());

        OpenMagneticsTesting::check_turns_description(coil);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Sections_Two_Windings_Together_One_Not) {
        settings->set_coil_wind_even_if_not_fit(true);
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

        CHECK_EQUAL(4, coil.get_sections_description()->size());
        CHECK_EQUAL(2, coil.get_sections_description().value()[0].get_partial_windings().size());
        CHECK_EQUAL(1, coil.get_sections_description().value()[2].get_partial_windings().size());
        CHECK_EQUAL("winding 0", coil.get_sections_description().value()[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL("winding 1", coil.get_sections_description().value()[0].get_partial_windings()[1].get_winding());
        CHECK_EQUAL("winding 2", coil.get_sections_description().value()[2].get_partial_windings()[0].get_winding());
        auto virtualFunctionalDescription = coil.virtualize_functional_description();
        CHECK_EQUAL(2, virtualFunctionalDescription.size());
        CHECK_EQUAL(numberTurns[0] + numberTurns[1], virtualFunctionalDescription[0].get_number_turns());
        CHECK_EQUAL(numberParallels[0], virtualFunctionalDescription[0].get_number_parallels());
        CHECK_EQUAL(numberTurns[2], virtualFunctionalDescription[1].get_number_turns());
        CHECK_EQUAL(numberParallels[2], virtualFunctionalDescription[1].get_number_parallels());

        OpenMagneticsTesting::check_turns_description(coil);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Layers_Two_Windings_Together) {
        settings->set_coil_wind_even_if_not_fit(true);
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
        CHECK_EQUAL(1, coil.get_layers_description()->size());
        CHECK_EQUAL(2, coil.get_layers_description().value()[0].get_partial_windings().size());
        CHECK_EQUAL("winding 0", coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL("winding 1", coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());

        OpenMagneticsTesting::check_turns_description(coil);
        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Layers_Two_Windings_Together_One_Not) {
        settings->set_coil_wind_even_if_not_fit(true);
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
        CHECK_EQUAL(4, coil.get_layers_description()->size());
        CHECK_EQUAL(2, coil.get_layers_description().value()[0].get_partial_windings().size());
        CHECK_EQUAL(1, coil.get_layers_description().value()[2].get_partial_windings().size());
        CHECK_EQUAL("winding 0", coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL("winding 1", coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());
        CHECK_EQUAL("winding 2", coil.get_layers_description().value()[2].get_partial_windings()[0].get_winding());

        OpenMagneticsTesting::check_turns_description(coil);
        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turns_Two_Windings_Together) {
        settings->set_coil_wind_even_if_not_fit(true);
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
        CHECK_EQUAL(1, coil.get_layers_description()->size());
        CHECK_EQUAL(2, coil.get_layers_description().value()[0].get_partial_windings().size());
        CHECK_EQUAL("winding 0", coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL("winding 1", coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());

        OpenMagneticsTesting::check_turns_description(coil);
        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }

    TEST(Test_Wind_By_Turns_Two_Windings_Together_One_Not) {
        settings->set_coil_wind_even_if_not_fit(true);
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
        CHECK_EQUAL(4, coil.get_layers_description()->size());
        CHECK_EQUAL(2, coil.get_layers_description().value()[0].get_partial_windings().size());
        CHECK_EQUAL("winding 0", coil.get_layers_description().value()[0].get_partial_windings()[0].get_winding());
        CHECK_EQUAL("winding 1", coil.get_layers_description().value()[0].get_partial_windings()[1].get_winding());
        CHECK_EQUAL("winding 2", coil.get_layers_description().value()[2].get_partial_windings()[0].get_winding());

        OpenMagneticsTesting::check_turns_description(coil);
        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
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
            settings->reset();
        }
    }
}