#include "WireWrapper.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

SUITE(Wire) {
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
    double max_error = 0.05;

    TEST(Sample_Wire) {
        auto wireFilePath = masPath + "samples/magnetic/wire/round/0.000140.json";
        std::ifstream json_file(wireFilePath);
        auto wireJson = json::parse(json_file);

        OpenMagnetics::WireWrapper wire(wireJson);

        auto conductingDiameter = wire.get_conducting_diameter().value().get_nominal().value();

        CHECK(conductingDiameter == wireJson["conductingDiameter"]["nominal"]);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(5.4e-05);
        double expectedValue = 0.755;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(1.1e-05);
        double expectedValue = 0.64;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00048);
        double expectedValue = 0.87;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(5.4e-05, 2);
        double expectedValue = 0.616;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(1.1e-05, 2);
        double expectedValue = 0.455;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00048, 2);
        double expectedValue = 0.8;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(5.4e-05, 3);
        double expectedValue = 0.523;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(1.1e-05, 3);
        double expectedValue = 0.334;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00048, 3);
        double expectedValue = 0.741;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(5.4e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.79;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(1.3e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.71;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00048, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.89;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(5.4e-05, 2, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.65;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(1.3e-05, 2, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.52;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00048, 2, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.81;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(5.4e-05, 3, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.55;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(4e-05, 3, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.51;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00048, 3, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.74;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Thin_Round_Insulated_Wire_1_Layer_Thin_Layer_Thickness) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.000101, 1, 3.81e-05, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.321961;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Thick_Round_Insulated_Wire_1_Layer_Thick_Layer_Thickness) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00129, 1, 7.62e-05, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.799184;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Outer_Diameter_Thick_Round_Insulated_Wire_1_Layer_Thick_Layer_Thickness) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_round(0.00129, 1, 7.62e-05, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00144;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Filling_Factors_Thick_Round_Insulated_Wire_3_Layer_Thick_Layer_Thickness) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_round(0.00129, 3, 7.62e-05, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.5446;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Outer_Diameter_Thick_Round_Insulated_Wire_3_Layers_Thick_Layer_Thickness) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_round(0.00129, 3, 7.62e-05, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.001748;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Filling_Factor_Thick_Litz_Wire_Served_1_Layer_Few_Strands) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_served_litz(0.0001, 66, 1, 1, OpenMagnetics::WireStandard::IEC_60317, false);
        double expectedValue = 0.458122;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Outer_Diameter_Litz_Wire_Grade_1_Served_1_Layer_Few_Strands) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_served_litz(2e-05, 10, 1, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.000126;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Litz_Wire_Grade_2_Served_1_Layer_1_Few_Strands) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_served_litz(2e-05, 10, 2, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.000137;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Thick_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Few_Strands) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_insulated_litz(0.000102, 66, 3, 7.62e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00152908;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Thin_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Many_Strands) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_insulated_litz(5.1e-05, 825, 3, 7.62e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00253238;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Thin_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Many_Strands_Diameter_Not_In_Db) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_outer_diameter_insulated_litz(5.42e-05, 825, 3, 7.62e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00253238;

        CHECK(expectedOuterDiameter < outerDiameter);
    }

    TEST(Test_Filling_Factor_Litz_Wire_Grade_2_Served_1_Layer_1_Few_Strands) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_filling_factor_served_litz(2e-05, 10, 2, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.235;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Filling_Factor_Thick_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Few_Strands) {
        auto outerDiameter = OpenMagnetics::WireWrapper::get_filling_factor_insulated_litz(0.000102, 66, 3, 7.62e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.3449;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Width_Small_Rectangular_Grade_1) {
        auto outerWidth = OpenMagnetics::WireWrapper::get_outer_width_rectangular(0.002, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterWidth = 0.00206;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Width_Small_Rectangular_Grade_2) {
        auto outerWidth = OpenMagnetics::WireWrapper::get_outer_width_rectangular(0.002, 2, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterWidth = 0.00217;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Width_Large_Rectangular_Grade_1) {
        auto outerWidth = OpenMagnetics::WireWrapper::get_outer_width_rectangular(0.016, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterWidth = 0.01608;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Width_Large_Rectangular_Grade_2) {
        auto outerWidth = OpenMagnetics::WireWrapper::get_outer_width_rectangular(0.016, 2, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterWidth = 0.01614;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Height_Small_Rectangular_Grade_1) {
        auto outerHeight = OpenMagnetics::WireWrapper::get_outer_height_rectangular(0.0008, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00088;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Height_Small_Rectangular_Grade_2) {
        auto outerHeight = OpenMagnetics::WireWrapper::get_outer_height_rectangular(0.0008, 2, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00092;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Height_Large_Rectangular_Grade_1) {
        auto outerHeight = OpenMagnetics::WireWrapper::get_outer_height_rectangular(0.0045, 1, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00456;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Height_Large_Rectangular_Grade_2) {
        auto outerHeight = OpenMagnetics::WireWrapper::get_outer_height_rectangular(0.0045, 2, OpenMagnetics::WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00467;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Filling_Factor_Small_Rectangular_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_rectangular(0.002, 0.0008, 2, OpenMagnetics::WireStandard::IEC_60317);
        double expectedFillingFactor = 0.720267;

        CHECK_CLOSE(expectedFillingFactor, fillingFactor, max_error * expectedFillingFactor);
    }

    TEST(Test_Filling_Factor_Large_Rectangular_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor_rectangular(0.016, 0.0045, 2, OpenMagnetics::WireStandard::IEC_60317);
        double expectedFillingFactor = 0.948615;

        CHECK_CLOSE(expectedFillingFactor, fillingFactor, max_error * expectedFillingFactor);
    }

    TEST(Test_Conducting_Area_Small_Rectangular) {
        auto conductingArea = OpenMagnetics::WireWrapper::get_conducting_area_rectangular(0.002, 0.0008, OpenMagnetics::WireStandard::IEC_60317);
        double expectedConductingArea = 0.00000146;

        CHECK_CLOSE(expectedConductingArea, conductingArea, max_error * expectedConductingArea);
    }

    TEST(Test_Conducting_Area_Large_Rectangular) {
        auto conductingArea = OpenMagnetics::WireWrapper::get_conducting_area_rectangular(0.016, 0.0045, OpenMagnetics::WireStandard::IEC_60317);
        double expectedConductingArea = 0.00007114;

        CHECK_CLOSE(expectedConductingArea, conductingArea, max_error * expectedConductingArea);
    }

}