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


SUITE(CoilSectionsDescription) {

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced) {
        std::vector<int64_t> numberTurns = {41};
        std::vector<int64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Full_Turns) {
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Full_Parallels) {
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Full_Parallels_Multiwinding) {
        std::vector<int64_t> numberTurns = {2, 5};
        std::vector<int64_t> numberParallels = {7, 7};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel);
    }

    TEST(Wind_By_Section_Wind_By_Consecutive_Parallels_Not_Balanced_Vertical) {
        std::vector<int64_t> numberTurns = {41};
        std::vector<int64_t> numberParallels = {3};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_0) {
        std::vector<int64_t> numberTurns = {9};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_1) {
        std::vector<int64_t> numberTurns = {6};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_2) {
        std::vector<int64_t> numberTurns = {5};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_3) {
        std::vector<int64_t> numberTurns = {5};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_4) {
        std::vector<int64_t> numberTurns = {91};
        std::vector<int64_t> numberParallels = {2};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 3;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_5) {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 7;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random_6) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {43};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 5;

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);

        OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, OpenMagnetics::WindingOrientation::VERTICAL);
    }

    TEST(Wind_By_Section_Random) {
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
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;

            auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
    }

    TEST(Wind_By_Section_Random_Multiwinding) {
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
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;
            if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                bobbinWidth *= numberTurns.size();
            }
            else {
                bobbinHeight *= numberTurns.size();
            }

            auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            OpenMagneticsTesting::check_sections_description(coil, numberTurns, numberParallels, interleavingLevel, windingOrientation);
        }
    }

    TEST(Wind_By_Section_With_Insulation_Sections) {
        std::vector<int64_t> numberTurns = {23, 42};
        std::vector<int64_t> numberParallels = {2, 1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;

        auto wires = std::vector<OpenMagnetics::WireWrapper>({OpenMagnetics::find_wire_by_name("0.014 - Grade 1")});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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

    TEST(Wind_By_Section_Pattern) {
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

    TEST(Wind_By_Layer_Wind_One_Section_One_Layer) {
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

    TEST(Wind_By_Layer_Wind_One_Section_Two_Layers) {
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

    TEST(Wind_By_Layer_Wind_One_Section_One_Layer_Two_Parallels) {
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

    TEST(Wind_By_Layer_Wind_One_Section_Two_Layers_Two_Parallels) {
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

    TEST(Wind_By_Layer_Wind_Two_Sections_Two_Layers_Two_Parallels) {
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

    TEST(Wind_By_Layer_Wind_Two_Sections_One_Layer_One_Parallel) {
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

    TEST(Wind_By_Layer_Wind_Two_Sections_One_Layer_Two_Parallels) {
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

    TEST(Wind_By_Layer_Wind_Two_Sections_Two_Layers_One_Parallel) {
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

    TEST(Wind_By_Layer_Wind_Vertical_Winding_Horizontal_Layers) {
        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Vertical_Winding_Vertical_Layers) {
        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Horizontal_Winding_Horizontal_Layers) {
        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Horizontal_Winding_Vertical_Layers) {
        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation, layersOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil, layersOrientation);
    }

    TEST(Wind_By_Layer_Wind_Horizontal_Winding) {
        std::vector<int64_t> numberTurns = {7};
        std::vector<int64_t> numberParallels = {1};
        double wireDiameter = 0.000509;
        int64_t numberMaximumTurnsPerLayer = 2;
        int64_t numberMaximumLayers = 2;
        uint8_t interleavingLevel = 1;
        double bobbinHeight = double(numberMaximumTurnsPerLayer) * wireDiameter;
        double bobbinWidth = double(numberMaximumLayers) * double(interleavingLevel) * 0.000509;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0}; 

        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
     
        OpenMagneticsTesting::check_layers_description(coil);
    }

    TEST(Wind_By_Layer_Random_0) {
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

    TEST(Wind_By_Layer_Random) {
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

    TEST(Wind_By_Layer_With_Insulation_Layers) {
        std::vector<int64_t> numberTurns = {23, 42};
        std::vector<int64_t> numberParallels = {2, 1};
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;

        auto wires = std::vector<OpenMagnetics::WireWrapper>({OpenMagnetics::find_wire_by_name("0.014 - Grade 1")});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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

    TEST(Wind_By_Turn_Wind_One_Section_One_Layer) {
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
    }

    TEST(Wind_By_Turn_Random_Multiwinding) {
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
            uint8_t interleavingLevel = std::rand() % 10 + 1;
            interleavingLevel = std::min(std::max(uint8_t(1U), uint8_t(numberPhysicalTurns)), interleavingLevel);
            auto windingOrientation = std::rand() % 2? OpenMagnetics::WindingOrientation::VERTICAL : OpenMagnetics::WindingOrientation::HORIZONTAL;
            if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
                bobbinWidth *= numberTurns.size();
            }
            else {
                bobbinHeight *= numberTurns.size();
            }

            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
            OpenMagneticsTesting::check_turns_description(coil);
        }
    }

    TEST(Wind_By_Turn_Random_Multiwinding_0) {
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
        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Wind_By_Turn_Random_Multiwinding_1) {
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
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        if (windingOrientation == OpenMagnetics::WindingOrientation::HORIZONTAL) {
            bobbinWidth *= numberTurns.size();
        }
        else {
            bobbinHeight *= numberTurns.size();
        }

        auto coil = OpenMagneticsTesting::get_quick_coil_no_compact(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);

        OpenMagneticsTesting::check_layers_description(coil);
        OpenMagneticsTesting::check_turns_description(coil);
    }

    TEST(Wind_By_Turn_Wind_One_Section_One_Layer_Rectangular_No_Bobbin) {
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

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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
    }

}
