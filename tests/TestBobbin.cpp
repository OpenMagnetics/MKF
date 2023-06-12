#include "BobbinWrapper.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <vector>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;
#include <typeinfo>


SUITE(Bobbin) {
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
    double max_error = 0.01;

    TEST(Sample_Bobbin) {
        auto wireFilePath = masPath + "samples/magnetic/bobbin/bobbin_E19_5.json";
        std::ifstream json_file(wireFilePath);
        auto bobbinJson = json::parse(json_file);

        OpenMagnetics::BobbinWrapper bobbin(bobbinJson);

        double expectedColumnThickness = 0.00080;
        double expectedWallThickness = 0.00080;

        CHECK_CLOSE(expectedColumnThickness, bobbin.get_processed_description().value().get_column_thickness(), max_error * expectedColumnThickness);
        CHECK_CLOSE(expectedWallThickness, bobbin.get_processed_description().value().get_wall_thickness(), max_error * expectedWallThickness);
    }

    TEST(Get_Filling_Factors_Bobbin_Medium) {
        auto fillingFactor = OpenMagnetics::BobbinWrapper::get_filling_factor(0.009, 0.0275);

        double expectedValue = 0.737;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Small) {
        auto fillingFactor = OpenMagnetics::BobbinWrapper::get_filling_factor(0.002, 0.005);

        double expectedValue = 0.515;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Large) {
        auto fillingFactor = OpenMagnetics::BobbinWrapper::get_filling_factor(0.019, 0.057);

        double expectedValue = 0.738;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Outside_Above) {
        auto fillingFactor = OpenMagnetics::BobbinWrapper::get_filling_factor(1, 1);

        double expectedValue = 0.738;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Bobbin_Outside_Below) {
        auto fillingFactor = OpenMagnetics::BobbinWrapper::get_filling_factor(0, 0);

        double expectedValue = 0.377;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Winding_Window_Dimensions_Medium) {
        auto windingWindowDimensions = OpenMagnetics::BobbinWrapper::get_winding_window_dimensions(0.012, 0.027);

        double expectedWidthValue = 0.00985;
        double expectedHeightValue = 0.0237;
        double width = windingWindowDimensions.first;
        double height = windingWindowDimensions.second;

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }

    TEST(Get_Winding_Window_Dimensions_Too_Small) {
        auto windingWindowDimensions = OpenMagnetics::BobbinWrapper::get_winding_window_dimensions(0.002, 0.002);

        double expectedWidthValue = 0.00116;
        double expectedHeightValue = 0.001;
        double width = windingWindowDimensions.first;
        double height = windingWindowDimensions.second;

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }

    TEST(Get_Winding_Window_Dimensions_Too_Large) {
        auto windingWindowDimensions = OpenMagnetics::BobbinWrapper::get_winding_window_dimensions(0.1, 0.1);

        double expectedWidthValue = 0.09575;
        double expectedHeightValue = 0.0943;
        double width = windingWindowDimensions.first;
        double height = windingWindowDimensions.second;

        CHECK_CLOSE(expectedWidthValue, width, max_error * expectedWidthValue);
        CHECK_CLOSE(expectedHeightValue, height, max_error * expectedHeightValue);
    }
}