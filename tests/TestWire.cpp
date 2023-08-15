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
    double max_error = 0.01;

    TEST(Sample_Wire) {
        auto wireFilePath = masPath + "samples/magnetic/wire/solid/0.000140.json";
        std::ifstream json_file(wireFilePath);
        auto wireJson = json::parse(json_file);

        OpenMagnetics::WireWrapper wire(wireJson);

        auto conductingDiameter = wire.get_conducting_diameter().value().get_nominal().value();

        CHECK(conductingDiameter == wireJson["conductingDiameter"]["nominal"]);
    }

    TEST(Get_Filling_Factors_Medium_Solid_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(5.4e-05);
        double expectedValue = 0.755;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Small_Solid_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(1.1e-05);
        double expectedValue = 0.64;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Large_Solid_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(0.00048);
        double expectedValue = 0.87;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Medium_Solid_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(5.4e-05, 2);
        double expectedValue = 0.616;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Small_Solid_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(1.1e-05, 2);
        double expectedValue = 0.455;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Large_Solid_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(0.00048, 2);
        double expectedValue = 0.8;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Medium_Solid_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(5.4e-05, 3);
        double expectedValue = 0.523;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Small_Solid_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(1.1e-05, 3);
        double expectedValue = 0.334;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Large_Solid_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(0.00048, 3);
        double expectedValue = 0.741;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Medium_Solid_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(5.4e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.79;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Small_Solid_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(1.3e-05, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.71;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Large_Solid_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(0.00048, 1, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.89;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Medium_Solid_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(5.4e-05, 2, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.65;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Small_Solid_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(1.3e-05, 2, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.52;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Large_Solid_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(0.00048, 2, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.81;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Medium_Solid_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(5.4e-05, 3, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.55;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Small_Solid_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(4e-05, 3, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.51;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Get_Filling_Factors_Large_Solid_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::WireWrapper::get_filling_factor(0.00048, 3, OpenMagnetics::WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.74;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }
}