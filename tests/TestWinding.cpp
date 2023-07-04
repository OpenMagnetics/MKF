#include "WindingWrapper.h"
#include "json.hpp"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <limits>
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
        size_t numberInsulationSections = 0;

        // json mierda;
        // to_json(mierda, windingWindow);
        for (auto& section : sectionsDescription){
            if(section.get_type() == OpenMagnetics::ElectricalType::INSULATION) {
                numberInsulationSections++;
                sectionsArea += section.get_dimensions()[0] * section.get_dimensions()[1];
            }
            else {
                // to_json(mierda, section);
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
                CHECK(section.get_filling_factor().value() > 0);
            }
        }
        for (size_t i = 0; i < sectionsDescription.size() - 1; ++i){
            if(sectionsDescription[i].get_type() == OpenMagnetics::ElectricalType::INSULATION) {
            }
            else {
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
        }

        CHECK(OpenMagnetics::roundFloat(bobbinArea, 6) == OpenMagnetics::roundFloat(sectionsArea, 6));
        for (size_t i = 0; i < numberAssignedParallels.size() - 1; ++i){
            CHECK(round(numberAssignedParallels[i]) == round(numberParallels[i]));
            CHECK(round(numberAssignedPhysicalTurns[i]) == round(numberTurns[i] * numberParallels[i]));
        }
        CHECK(sectionsDescription.size() - numberInsulationSections == (interleavingLevel * numberTurns.size()));
        CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced) {
        std::vector<uint64_t> numberTurns = {41};
        std::vector<uint64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Full_Turns) {
        std::vector<uint64_t> numberTurns = {2};
        std::vector<uint64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Full_Parallels) {
        std::vector<uint64_t> numberTurns = {2};
        std::vector<uint64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 7;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Full_Parallels_Multiwinding) {
        std::vector<uint64_t> numberTurns = {2, 5};
        std::vector<uint64_t> numberParallels = {7, 7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 7;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel);

    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced_Vertical) {
        std::vector<uint64_t> numberTurns = {41};
        std::vector<uint64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 2;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

    }

    TEST(Wind_By_Section_Random_0) {
        std::vector<uint64_t> numberTurns = {9};
        std::vector<uint64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_1) {
        std::vector<uint64_t> numberTurns = {6};
        std::vector<uint64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_2) {
        std::vector<uint64_t> numberTurns = {5};
        std::vector<uint64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_3) {
        std::vector<uint64_t> numberTurns = {5};
        std::vector<uint64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_4) {
        std::vector<uint64_t> numberTurns = {91};
        std::vector<uint64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 3;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_5) {
        std::vector<uint64_t> numberTurns = {23};
        std::vector<uint64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 7;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_6) {
        std::vector<uint64_t> numberTurns = {1};
        std::vector<uint64_t> numberParallels = {43};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 5;

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

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
            uint64_t numberPhysicalTurns = numberTurns[0] * numberParallels[0];
            uint64_t interleavingLevel = std::rand() % 10 + 1;
            interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;

            auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
    }

    TEST(Wind_By_Section_Random_Multiwinding) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<uint64_t> numberTurns;
            std::vector<uint64_t> numberParallels;
            uint64_t numberPhysicalTurns = std::numeric_limits<uint64_t>::max();
            for (size_t windingIndex = 0; windingIndex < std::rand() % 10 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1UL);
                numberParallels.push_back(std::rand() % 100 + 1UL);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
            }
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            uint64_t interleavingLevel = std::rand() % 10 + 1;
            interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;
            if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                bobbinWidth *= numberTurns.size();
            }
            else {
                bobbinHeight *= numberTurns.size();
            }

            auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
    }

    TEST(Wind_By_Section_With_Insulation_Sections) {
        std::vector<uint64_t> numberTurns = {23, 42};
        std::vector<uint64_t> numberParallels = {2, 1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint64_t interleavingLevel = 2;

        auto wires = std::vector<OpenMagnetics::WireWrapper>({OpenMagnetics::find_wire_by_name("0.014 - Grade 1")});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        double voltagePeakToPeak = 400;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        winding.set_inputs(inputs);
        winding.wind();
        auto log = winding.read_log();

        quick_check_sections_description(winding, numberTurns, numberParallels, interleavingLevel, sectionOrientation);
    }
}

SUITE(WindingLayersDescription) {

    void quick_check_layers_description(OpenMagnetics::WindingWrapper winding,
                                          OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL) {
        if (!winding.get_layers_description()) {
            return;
        }
        auto sections = winding.get_sections_description().value();
        std::map<std::string, std::vector<double>> dimensionsByName;
        std::map<std::string, std::vector<double>> coordinatesByName;

        for (auto& section : sections){
            auto layers = winding.get_layers_by_section(section.get_name());
            if(section.get_type() == OpenMagnetics::ElectricalType::INSULATION) {
            }
            else {
                auto sectionParallelsProportionExpected = section.get_partial_windings()[0].get_parallels_proportion();
                std::vector<double> sectionParallelsProportion(sectionParallelsProportionExpected.size(), 0);
                for (auto& layer : layers){
                    for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                        sectionParallelsProportion[i] += layer.get_partial_windings()[0].get_parallels_proportion()[i];
                    }
                    CHECK(layer.get_filling_factor().value() > 0);

                    dimensionsByName[layer.get_name()] = layer.get_dimensions();
                    coordinatesByName[layer.get_name()] = layer.get_coordinates();
                }
                for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                    CHECK(OpenMagnetics::roundFloat(sectionParallelsProportion[i], 9) == OpenMagnetics::roundFloat(sectionParallelsProportionExpected[i], 9));
                }
                for (size_t i = 0; i < layers.size() - 1; ++i){
                    if (layersOrientation == OpenMagnetics::WindingOrientation::VERTICAL) {
                        CHECK(layers[i].get_coordinates()[0] < layers[i + 1].get_coordinates()[0]);
                        CHECK(layers[i].get_coordinates()[1] == layers[i + 1].get_coordinates()[1]);
                        CHECK(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                    } 
                    else if (layersOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                        CHECK(layers[i].get_coordinates()[1] > layers[i + 1].get_coordinates()[1]);
                        CHECK(layers[i].get_coordinates()[0] == layers[i + 1].get_coordinates()[0]);
                        CHECK(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                    }
                }
            }

        }

        CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
    }


    TEST(Wind_By_Layer_Wind_One_Section_One_Layer) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 9;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Wind_One_Section_Two_Layers) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 6;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }
    TEST(Wind_By_Layer_Wind_One_Section_One_Layer_Two_Parallels) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {2};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 15;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Wind_One_Section_Two_Layers_Two_Parallels) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {2};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 6;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }
    TEST(Wind_By_Layer_Wind_Two_Sections_Two_Layers_Two_Parallels) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {2};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 6;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 2;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Wind_Two_Sections_One_Layer_One_Parallel) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 6;
        uint64_t numberMaximumLayers = 1;
        uint64_t interleavingLevel = 2;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Wind_Two_Sections_One_Layer_Two_Parallels) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {2};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 6;
        uint64_t numberMaximumLayers = 1;
        uint64_t interleavingLevel = 2;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        auto layersDescription = winding.get_layers_description().value();
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Wind_Two_Sections_Two_Layers_One_Parallel) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 2;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 2;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
     
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Wind_Vertical_Winding_Horizontal_Layers) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 2;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        quick_check_layers_description(winding, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Vertical_Winding_Vertical_Layers) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 2;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        quick_check_layers_description(winding, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Horizontal_Winding_Horizontal_Layers) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 2;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        quick_check_layers_description(winding, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Horizontal_Winding_Vertical_Layers) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 2;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        quick_check_layers_description(winding, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Horizontal_Winding) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 2;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
     
        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Random_0) {
        std::vector<uint64_t> numberTurns = {5};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 1;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 2;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        quick_check_layers_description(winding);
    }

    TEST(Wind_By_Layer_Random) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<uint64_t> numberTurns = {std::rand() % 10 + 1UL};
            std::vector<uint64_t> numberParallels = {std::rand() % 3 + 1U};
            double wireDiameter = 0.000509;
            uint64_t numberMaximumTurnsPerLayer = std::rand() % 4 + 1U;
            uint64_t numberMaximumLayers = std::rand() % 3 + 1U;
            uint64_t interleavingLevel = std::rand() % 10 + 1;
            uint64_t numberPhysicalTurns = numberTurns[0] * numberParallels[0];
            interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
            double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
            double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            // std::cout << "numberTurns: " << numberTurns[0] << std::endl;
            // std::cout << "numberParallels: " << numberParallels[0] << std::endl;
            // std::cout << "numberMaximumTurnsPerLayer: " << numberMaximumTurnsPerLayer << std::endl;
            // std::cout << "numberMaximumLayers: " << numberMaximumLayers << std::endl;
            // std::cout << "interleavingLevel: " << interleavingLevel << std::endl;

            auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
            quick_check_layers_description(winding);
        }
    }

    TEST(Wind_By_Layer_With_Insulation_Layers) {
        std::vector<uint64_t> numberTurns = {23, 42};
        std::vector<uint64_t> numberParallels = {2, 1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint64_t interleavingLevel = 2;

        auto wires = std::vector<OpenMagnetics::WireWrapper>({OpenMagnetics::find_wire_by_name("0.014 - Grade 1")});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        double voltagePeakToPeak = 400;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        winding.set_inputs(inputs);
        winding.wind();
        auto log = winding.read_log();

        quick_check_layers_description(winding);
    }
}

SUITE(WindingTurnsDescription) {

    void quick_check_layers_description(OpenMagnetics::WindingWrapper winding,
                                          OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL) {
        if (!winding.get_layers_description()) {
            return;
        }
        auto sections = winding.get_sections_description().value();
        std::map<std::string, std::vector<double>> dimensionsByName;
        std::map<std::string, std::vector<double>> coordinatesByName;

        for (auto& section : sections){
            auto layers = winding.get_layers_by_section(section.get_name());
            if(section.get_type() == OpenMagnetics::ElectricalType::INSULATION) {
            }
            else {
                auto sectionParallelsProportionExpected = section.get_partial_windings()[0].get_parallels_proportion();
                std::vector<double> sectionParallelsProportion(sectionParallelsProportionExpected.size(), 0);
                for (auto& layer : layers){
                    for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                        sectionParallelsProportion[i] += layer.get_partial_windings()[0].get_parallels_proportion()[i];
                    }
                    CHECK(layer.get_filling_factor().value() > 0);

                    dimensionsByName[layer.get_name()] = layer.get_dimensions();
                    coordinatesByName[layer.get_name()] = layer.get_coordinates();
                }
                for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                    CHECK(OpenMagnetics::roundFloat(sectionParallelsProportion[i], 9) == OpenMagnetics::roundFloat(sectionParallelsProportionExpected[i], 9));
                }
                for (size_t i = 0; i < layers.size() - 1; ++i){
                    if (layersOrientation == OpenMagnetics::WindingOrientation::VERTICAL) {
                        CHECK(layers[i].get_coordinates()[0] < layers[i + 1].get_coordinates()[0]);
                        CHECK(layers[i].get_coordinates()[1] == layers[i + 1].get_coordinates()[1]);
                        CHECK(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                    } 
                    else if (layersOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                        CHECK(layers[i].get_coordinates()[1] > layers[i + 1].get_coordinates()[1]);
                        CHECK(layers[i].get_coordinates()[0] == layers[i + 1].get_coordinates()[0]);
                        CHECK(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                    }
                }
            }

        }

        CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
    }


    void quick_check_turns_description(OpenMagnetics::WindingWrapper winding) {
        if (!winding.get_turns_description()) {
            return;
        }

        std::vector<std::vector<double>> parallelProportion;
        for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
            parallelProportion.push_back(std::vector<double>(winding.get_number_parallels(windingIndex), 0));
        }

        auto turns = winding.get_turns_description().value();
        std::map<std::string, std::vector<double>> dimensionsByName;
        std::map<std::string, std::vector<double>> coordinatesByName;

        for (auto& turn : turns){
            auto windingIndex = winding.get_winding_index_by_name(turn.get_winding());
            parallelProportion[windingIndex][turn.get_parallel()] += 1.0 / winding.get_number_turns(windingIndex);
            dimensionsByName[turn.get_name()] = turn.get_dimensions().value();
            coordinatesByName[turn.get_name()] = turn.get_coordinates();
        // json mierda;
        // to_json(mierda, turn);

        // std::cout << "mierda: " << mierda << std::endl;

        }

        for (size_t windingIndex = 0; windingIndex < winding.get_functional_description().size(); ++windingIndex) {
            for (size_t parallelIndex = 0; parallelIndex < winding.get_number_parallels(windingIndex); ++parallelIndex) {
                CHECK(OpenMagnetics::roundFloat(parallelProportion[windingIndex][parallelIndex], 9) == 1);
            }
        }
        CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
    }


    TEST(Wind_By_Turn_Wind_One_Section_One_Layer) {
        std::vector<uint64_t> numberTurns = {7};
        std::vector<uint64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        uint64_t numberMaximumTurnsPerLayer = 9;
        uint64_t numberMaximumLayers = 2;
        uint64_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);
        quick_check_turns_description(winding);
    }

    TEST(Wind_By_Turn_Random_Multiwinding) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            std::vector<uint64_t> numberTurns;
            std::vector<uint64_t> numberParallels;
            uint64_t numberPhysicalTurns = std::numeric_limits<uint64_t>::max();
            for (size_t windingIndex = 0; windingIndex < std::rand() % 10 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1UL);
                numberParallels.push_back(std::rand() % 100 + 1UL);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
            }
            double bobbinHeight = 0.01;
            double bobbinWidth = 0.01;
            std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
            uint64_t interleavingLevel = std::rand() % 10 + 1;
            interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;
            if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                bobbinWidth *= numberTurns.size();
            }
            else {
                bobbinHeight *= numberTurns.size();
            }

            auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            quick_check_turns_description(winding);
        }
    }

    TEST(Wind_By_Turn_Random_Multiwinding_0) {
        std::vector<uint64_t> numberTurns;
        std::vector<uint64_t> numberParallels;
        uint64_t numberPhysicalTurns = std::numeric_limits<uint64_t>::max();
        for (size_t windingIndex = 0; windingIndex < 1UL; ++windingIndex)
        {
            numberTurns.push_back(4);
            numberParallels.push_back(12);
            numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
        }
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 10;
        interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        quick_check_turns_description(winding);
    }

    TEST(Wind_By_Turn_Random_Multiwinding_1) {
        std::vector<uint64_t> numberTurns = {80};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t numberPhysicalTurns = std::numeric_limits<uint64_t>::max();

        for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
        {
            numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
        }
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint64_t interleavingLevel = 9;
        interleavingLevel = std::min(numberPhysicalTurns, interleavingLevel);
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto winding = OpenMagneticsTesting::get_quick_winding_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        quick_check_layers_description(winding);
        quick_check_turns_description(winding);
    }

}