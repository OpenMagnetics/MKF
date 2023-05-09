#include "WindingWrapper.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <vector>
#include <time.h>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;
#include <typeinfo>

SUITE(WindingFunctionalDescription) {
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");

    TEST(Inductor_42_Turns) {
        auto windingFilePath = masPath + "samples/magnetic/winding/inductor_42_turns.json";
        std::ifstream json_file(windingFilePath);

        auto windingJson = json::parse(json_file);

        OpenMagnetics::WindingWrapper winding(windingJson);

        auto function_description = winding.get_functional_description()[0];

        json windingWrapperJson;
        to_json(windingWrapperJson, function_description);

        CHECK(windingWrapperJson == windingJson["functionalDescription"][0]);
    }
}

SUITE(WindingSectionsDescription) {

    OpenMagnetics::WindingWrapper get_quick_winding(std::vector<uint64_t> numberTurns,
                                                    std::vector<uint64_t> numberParallels,
                                                    double bobbinHeight,
                                                    double bobbinWidth,
                                                    std::vector<double> bobbinCenterCoodinates,
                                                    uint64_t interleavingLevel,
                                                    OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL){
        json windingJson;
        windingJson["functionalDescription"] = json::array();

        windingJson["bobbin"] = json();
        windingJson["bobbin"]["processedDescription"] = json();
        windingJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
        windingJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
        windingJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

        json windingWindow;
        windingWindow["height"] = bobbinHeight;
        windingWindow["width"] = bobbinWidth;
        windingWindow["coordinates"] = json::array();
        for (auto& coord : bobbinCenterCoodinates) {
            windingWindow["coordinates"].push_back(coord);
        }
        windingJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

        for (size_t i = 0; i < numberTurns.size(); ++i){
            json individualWindingJson;
            individualWindingJson["name"] = "winding " + std::to_string(i);
            individualWindingJson["numberTurns"] = numberTurns[i];
            individualWindingJson["numberParallels"] = numberParallels[i];
            individualWindingJson["isolationSide"] = "primary";
            individualWindingJson["wire"] = "Litz 450x01";
            windingJson["functionalDescription"].push_back(individualWindingJson);
        }

        OpenMagnetics::WindingWrapper winding(windingJson, interleavingLevel, windingOrientation);
        return winding;
    }

    void quick_check_sections_description(OpenMagnetics::WindingWrapper winding,
                                          std::vector<uint64_t> numberTurns,
                                          std::vector<uint64_t> numberParallels,
                                          uint64_t interleavingLevel,
                                          OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL) {

        auto bobbin = std::get<OpenMagnetics::Bobbin>(winding.get_bobbin());
        auto windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
        auto bobbinArea = windingWindow.get_width().value() * windingWindow.get_height().value();
        auto sectionsDescription = winding.get_sections_description().value();
        std::vector<double> numberAssignedParallels(numberTurns.size(), 0);
        std::vector<uint64_t> numberAssignedPhysicalTurns(numberTurns.size(), 0);
        std::map<std::string, std::vector<double>> dimensionsByName;
        std::map<std::string, std::vector<double>> coordinatesByName;
        double sectionsArea = 0;

        // json mierda;
        // to_json(mierda, windingWindow);
        // std::cout << mierda << std::endl;
        for (auto& section : sectionsDescription){
            // to_json(mierda, section);
            // std::cout << mierda << std::endl;
            for (auto& partialWinding : section.get_partial_windings()) {
                auto currentIndividualWindingIndex = winding.get_winding_index_by_name(partialWinding.get_winding());
                auto currentIndividualWinding = winding.get_winding_by_name(partialWinding.get_winding());
                sectionsArea += section.get_dimensions()[0] * section.get_dimensions()[1];
                dimensionsByName[section.get_name()] = section.get_dimensions();
                coordinatesByName[section.get_name()] = section.get_coordinates();
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, 6) >= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[0] - windingWindow.get_width().value() / 2, 6));
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[0] + section.get_dimensions()[0] / 2, 6) <= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2, 6));
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[1] - section.get_dimensions()[1] / 2, 6) >= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[1] - windingWindow.get_height().value() / 2, 6));
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[1] + section.get_dimensions()[1] / 2, 6) <= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[1] + windingWindow.get_height().value() / 2, 6));

                for (auto& parallelProportion : partialWinding.get_parallels_proportion()) {
                    numberAssignedParallels[currentIndividualWindingIndex] += parallelProportion;
                    numberAssignedPhysicalTurns[currentIndividualWindingIndex] += round(parallelProportion * currentIndividualWinding.get_number_turns());
                }
            }
        }
        for (size_t i = 0; i < sectionsDescription.size() - 1; ++i){
            if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                CHECK(sectionsDescription[i].get_coordinates()[0] < sectionsDescription[i + 1].get_coordinates()[0]);
                CHECK(sectionsDescription[i].get_coordinates()[1] == sectionsDescription[i + 1].get_coordinates()[1]);
                CHECK(sectionsDescription[i].get_coordinates()[2] == sectionsDescription[i + 1].get_coordinates()[2]);
            } 
            else if (windingOrientation == OpenMagnetics::WindingOrientation::VERTICAL) {
                CHECK(sectionsDescription[i].get_coordinates()[1] > sectionsDescription[i + 1].get_coordinates()[1]);
                CHECK(sectionsDescription[i].get_coordinates()[0] == sectionsDescription[i + 1].get_coordinates()[0]);
                CHECK(sectionsDescription[i].get_coordinates()[2] == sectionsDescription[i + 1].get_coordinates()[2]);
            }
        }

        // std::cout << "numberAssignedPhysicalTurns" << std::endl;
        // std::cout << numberAssignedPhysicalTurns << std::endl;
        // std::cout << "round(numberTurns * numberParallels)" << std::endl;
        // std::cout << round(numberTurns * numberParallels) << std::endl;
        CHECK(OpenMagnetics::roundFloat(bobbinArea, 6) == OpenMagnetics::roundFloat(sectionsArea, 6));
        for (size_t i = 0; i < numberAssignedParallels.size() - 1; ++i){
            CHECK(round(numberAssignedParallels[i]) == round(numberParallels[i]));
            CHECK(round(numberAssignedPhysicalTurns[i]) == round(numberTurns[i] * numberParallels[i]));
        }
        CHECK(sectionsDescription.size() == (interleavingLevel * numberTurns.size()));
        CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced) {
        std::vector<uint64_t> numberTurns = {41};
        std::vector<uint64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Full_Turns) {
        std::vector<uint64_t> numberTurns = {2};
        std::vector<uint64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Full_Parallels) {
        std::vector<uint64_t> numberTurns = {2};
        std::vector<uint64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 7;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Full_Parallels_Multiwinding) {
        std::vector<uint64_t> numberTurns = {2, 5};
        std::vector<uint64_t> numberParallels = {7, 7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 7;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced_Vertical) {
        std::vector<uint64_t> numberTurns = {41};
        std::vector<uint64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

    }

    TEST(Wind_By_Section_Random_0) {
        std::vector<uint64_t> numberTurns = {9};
        std::vector<uint64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_1) {
        std::vector<uint64_t> numberTurns = {6};
        std::vector<uint64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_2) {
        std::vector<uint64_t> numberTurns = {5};
        std::vector<uint64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_3) {
        std::vector<uint64_t> numberTurns = {5};
        std::vector<uint64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_4) {
        std::vector<uint64_t> numberTurns = {91};
        std::vector<uint64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<uint64_t> numberTurns = {std::rand() % 100 + 1UL};
            std::vector<uint64_t> numberParallels = {std::rand() % 100 + 1UL};
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            uint64_t interleavingLevel = std::rand() % 10 + 1;
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;

            auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
    }

    TEST(Wind_By_Section_Random_Multiwinding) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<uint64_t> numberTurns;
            std::vector<uint64_t> numberParallels;
            for (size_t windingIndex = 0; windingIndex < std::rand() % 10 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1UL);
                numberParallels.push_back(std::rand() % 100 + 1UL);
            }
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            uint64_t interleavingLevel = std::rand() % 10 + 1;
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;

            auto winding = get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
    }

    // TEST(Wind_By_Section_Transformer) {
    //     json windingJson;
    //     windingJson["functionalDescription"] = json::array();
    //     json primaryWindingJson;
    //     primaryWindingJson["name"] = "primary";
    //     primaryWindingJson["numberTurns"] = 42;
    //     primaryWindingJson["numberParallels"] = 3;
    //     primaryWindingJson["isolationSide"] = "primary";
    //     primaryWindingJson["wire"] = "Litz 450x01";

    //     // json secondaryWindingJson;
    //     // secondaryWindingJson = json();
    //     // secondaryWindingJson["name"] = "secondary";
    //     // secondaryWindingJson["numberTurns"] = 23;
    //     // secondaryWindingJson["numberParallels"] = 5;
    //     // secondaryWindingJson["isolationSide"] = "secondary";
    //     // secondaryWindingJson["wire"] = "Litz 450x01";

    //     windingJson["bobbin"] = json();
    //     windingJson["bobbin"]["processedDescription"] = json();
    //     windingJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
    //     windingJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
    //     windingJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    //     json windingWindow;
    //     windingWindow["height"] = 0.01;
    //     windingWindow["width"] = 0.01;
    //     windingJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    //     windingJson["functionalDescription"].push_back(primaryWindingJson);
    //     // windingJson["functionalDescription"].push_back(secondaryWindingJson);

    //     std::cout << "windingJson[functionalDescription].size()" << std::endl;
    //     std::cout << windingJson["functionalDescription"].size() << std::endl;
    //     OpenMagnetics::WindingWrapper winding(windingJson, 2);

    //     auto function_description = winding.get_functional_description()[0];
    // }
}