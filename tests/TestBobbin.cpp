#include "constructive_models/Bobbin.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;


SUITE(Bobbin) {
    auto masPath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("MAS/").string();
    double max_error = 0.05;

    TEST(Sample_Bobbin) {
        auto wireFilePath = masPath + "samples/magnetic/bobbin/bobbin_E19_5.json";
        std::ifstream json_file(wireFilePath);
        auto bobbinJson = json::parse(json_file);

        OpenMagnetics::Bobbin bobbin(bobbinJson);

        double expectedColumnThickness = 0.00080;
        double expectedWallThickness = 0.00080;

        CHECK_CLOSE(expectedColumnThickness, bobbin.get_processed_description().value().get_column_thickness(), max_error * expectedColumnThickness);
        CHECK_CLOSE(expectedWallThickness, bobbin.get_processed_description().value().get_wall_thickness(), max_error * expectedWallThickness);
    }

    TEST(Get_Filling_Factors_Bobbin_Medium) {
        auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.009, 0.0275);

        double expectedValue = 0.715;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Small) {
        auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.002, 0.005);

        double expectedValue = 0.53;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Large) {
        auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.019, 0.057);

        double expectedValue = 0.725;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Outside_Above) {
        auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(1, 1);

        double expectedValue = 0.738;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Outside_Below) {
        auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0, 0);

        double expectedValue = 0.377;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Winding_Window_Dimensions_Medium) {
        auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.012, 0.027);

        double expectedWidthValue = 0.00985;
        double expectedHeightValue = 0.02335;
        double width = windingWindowDimensions[0];
        double height = windingWindowDimensions[1];

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }

    TEST(Get_Winding_Window_Dimensions_Too_Small) {
        auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.001, 0.002);

        double expectedWidthValue = 0.0005;
        double expectedHeightValue = 0.0011;
        double width = windingWindowDimensions[0];
        double height = windingWindowDimensions[1];

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }

    TEST(Get_Winding_Window_Dimensions_Too_Large) {
        auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.1, 0.1);

        double expectedWidthValue = 0.0951;
        double expectedHeightValue = 0.0943;
        double width = windingWindowDimensions[0];
        double height = windingWindowDimensions[1];

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }

    TEST(Get_Winding_Window_Dimensions_Error) {
        auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.003325, 0.0108);

        double expectedWidthValue =  0.00245;
        double expectedHeightValue = 0.0094;
        double width = windingWindowDimensions[0];
        double height = windingWindowDimensions[1];

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }

    TEST(Get_Winding_Window_Dimensions_E_51) {
        auto core = OpenMagneticsTesting::get_quick_core("ER 51/10/38", json::parse("[]"), 1, "Dummy");
        auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
        auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(coreWindingWindow.get_width().value(), coreWindingWindow.get_height().value());

        auto widthThickness = coreWindingWindow.get_width().value() - windingWindowDimensions[0];
        auto heightThickness = (coreWindingWindow.get_height().value() - windingWindowDimensions[1]) / 2;

        CHECK(heightThickness <= widthThickness * 1.2);
    }

    TEST(Get_Winding_Window_Dimensions_All_Shapes_With_Bobbin) {
        settings->set_use_toroidal_cores(true);
        auto shapeNames = get_shape_names();
        for (auto shapeName : shapeNames) {
            if (shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI ")) {
                continue;
            }
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
            auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
            auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(coreWindingWindow.get_width().value(), coreWindingWindow.get_height().value());
        }
    }
}