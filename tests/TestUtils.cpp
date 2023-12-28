#include "InitialPermeability.h"
#include "Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <cfloat>
#include <limits>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

SUITE(Utils) {
    TEST(LoadDatabaseJson) {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/masData.json");

        std::ifstream ifs(masPath);
        json masData = json::parse(ifs);

        OpenMagnetics::load_databases(masData, true);
    }

    TEST(Modified_Bessel) {
        double calculatedValue = OpenMagnetics::modified_bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 1.2660658777520084;
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }


    TEST(Bessel) {
        double calculatedValue = OpenMagnetics::bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 0.7651976865579666;
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);

        double calculatedBerValue = OpenMagnetics::kelvin_function_real(0.0, 1.0);
        double expectedBerValue = 0.98438178;
        CHECK_CLOSE(expectedBerValue, calculatedBerValue, expectedBerValue * 0.001);

        double calculatedBeiValue = OpenMagnetics::kelvin_function_imaginary(0.0, 1.0);
        double expectedBeiValue = 0.24956604;
        CHECK_CLOSE(expectedBeiValue, calculatedBeiValue, expectedBeiValue * 0.001);

        double calculatedBerpValue = OpenMagnetics::derivative_kelvin_function_real(0.0, 1.0);
        double expectedBerpValue = -0.06244575217903096;
        CHECK_CLOSE(expectedBerpValue, calculatedBerpValue, fabs(expectedBerpValue) * 0.001);

        double calculatedBeipValue = OpenMagnetics::derivative_kelvin_function_imaginary(0.0, 1.0);
        double expectedBeipValue = 0.49739651146809727;
        CHECK_CLOSE(expectedBeipValue, calculatedBeipValue, expectedBeipValue * 0.001);

    }

    TEST(Test_Complete_Ellipitical_1_0) {
        double calculatedValue = OpenMagnetics::comp_ellint_1(0);
        double expectedValue = std::comp_ellint_1(0);
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_1_1) {
        double calculatedValue = OpenMagnetics::comp_ellint_1(1);
        double expectedValue = std::comp_ellint_1(1);
        CHECK(std::isnan(calculatedValue));
        CHECK(std::isnan(expectedValue));
    }

    TEST(Test_Complete_Ellipitical_1_2) {
        double calculatedValue = OpenMagnetics::comp_ellint_1(std::sin(std::numbers::pi / 18 / 2));
        double expectedValue = std::comp_ellint_1(std::sin(std::numbers::pi / 18 / 2));
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_2_0) {
        double calculatedValue = OpenMagnetics::comp_ellint_2(0);
        double expectedValue = std::comp_ellint_2(0);
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_2_1) {
        double calculatedValue = OpenMagnetics::comp_ellint_2(1);
        double expectedValue = std::comp_ellint_2(1);
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_2_2) {
        double calculatedValue = OpenMagnetics::comp_ellint_2(std::sin(std::numbers::pi / 18 / 2));
        double expectedValue = std::comp_ellint_2(std::sin(std::numbers::pi / 18 / 2));
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }
}