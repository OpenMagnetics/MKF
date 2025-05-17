#include "physical_models/CoreTemperature.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

std::map<CoreTemperatureModels, double> maximumAdmittedErrorTemperature = {
    {CoreTemperatureModels::KAZIMIERCZUK, 0.6}, {CoreTemperatureModels::MANIKTALA, 0.6},
    {CoreTemperatureModels::TDK, 0.71},         {CoreTemperatureModels::DIXON, 0.59},
    {CoreTemperatureModels::AMIDON, 0.6},
};
std::map<CoreTemperatureModels, std::vector<double>> testCoreTemperatureAverageErrors = {};
std::map<CoreTemperatureModels, double> testCoreTemperatureMaximumErrors = {};


double run_test_core_temperature(const CoreTemperatureModels& modelName,
                                 const std::string& shapeName,
                                 const std::string& materialName,
                                 double coreLosses,
                                 double ambientTemperature,
                                 double expectedCoreTemperature) {
    double maximumAdmittedErrorTemperatureValue = maximumAdmittedErrorTemperature[modelName];
    Core core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, materialName);
    auto coreTemperatureModel = CoreTemperatureModel::factory(modelName);

    auto coreTemperature = coreTemperatureModel->get_core_temperature(core, coreLosses, ambientTemperature);
    double calculatedTemperature = coreTemperature.get_maximum_temperature();

    double error = abs(expectedCoreTemperature - calculatedTemperature) / expectedCoreTemperature;

    if (error > maximumAdmittedErrorTemperatureValue) {
        std::cout << "error " << error * 100 << " %" << std::endl;
    }
    if (testCoreTemperatureMaximumErrors[modelName] < error) {
        testCoreTemperatureMaximumErrors[modelName] = error;
    }
    CHECK_CLOSE(calculatedTemperature, expectedCoreTemperature,
                expectedCoreTemperature * maximumAdmittedErrorTemperatureValue);

    return error;
}

void test_core_temperature_sotiris_47(CoreTemperatureModels modelName) {
    double meanError = 0;

    std::string coreShape = "ETD 49";
    std::string coreMaterial = "3C97";
    double coreLosses = 1.44;
    double ambientTemperature = 25;
    double expectedTemperature = 59;

    meanError += run_test_core_temperature(modelName, coreShape, coreMaterial, coreLosses, ambientTemperature,
                                           expectedTemperature);

    testCoreTemperatureAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core Temperature for " << coreMaterial << " with Model "
                  << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testCoreTemperatureAverageErrors[modelName].begin(),
                                 testCoreTemperatureAverageErrors[modelName].end()) /
                         testCoreTemperatureAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testCoreTemperatureMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_temperature_sotiris_46(CoreTemperatureModels modelName) {
    double meanError = 0;

    std::string coreShape = "ETD 44";
    std::string coreMaterial = "3C97";
    double coreLosses = 2.7;
    double ambientTemperature = 25;
    double expectedTemperature = 79;

    meanError += run_test_core_temperature(modelName, coreShape, coreMaterial, coreLosses, ambientTemperature,
                                           expectedTemperature);

    testCoreTemperatureAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core Temperature for " << coreMaterial << " with Model "
                  << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testCoreTemperatureAverageErrors[modelName].begin(),
                                 testCoreTemperatureAverageErrors[modelName].end()) /
                         testCoreTemperatureAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testCoreTemperatureMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_temperature_sotiris_40(CoreTemperatureModels modelName) {
    double meanError = 0;

    std::string coreShape = "ETD 29";
    std::string coreMaterial = "3C97";
    double coreLosses = 0.4;
    double ambientTemperature = 25;
    double expectedTemperature = 25 + 57;

    meanError += run_test_core_temperature(modelName, coreShape, coreMaterial, coreLosses, ambientTemperature,
                                           expectedTemperature);

    testCoreTemperatureAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core Temperature for " << coreMaterial << " with Model "
                  << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testCoreTemperatureAverageErrors[modelName].begin(),
                                 testCoreTemperatureAverageErrors[modelName].end()) /
                         testCoreTemperatureAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testCoreTemperatureMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_temperature_sotiris_37(CoreTemperatureModels modelName) {
    double meanError = 0;

    std::string coreShape = "EQ 25/6";
    std::string coreMaterial = "3C95";
    double ambientTemperature = 25;
    std::vector<double> coreLosses = {0.53, 0.76, 1.14, 1.49, 0.61, 0.88, 1.26, 1.58};
    std::vector<double> expectedTemperature = {44.5, 53.1, 67.3, 88.2, 47.1, 55.2, 67.6, 79.2};

    for (size_t i = 0; i < coreLosses.size(); ++i) {
        meanError += run_test_core_temperature(modelName, coreShape, coreMaterial, coreLosses[i], ambientTemperature,
                                               expectedTemperature[i]);
    }

    meanError /= coreLosses.size();
    testCoreTemperatureAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core Temperature for " << coreMaterial << " with Model "
                  << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testCoreTemperatureAverageErrors[modelName].begin(),
                                 testCoreTemperatureAverageErrors[modelName].end()) /
                         testCoreTemperatureAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testCoreTemperatureMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_temperature_miserable_40(CoreTemperatureModels modelName) {
    double meanError = 0;

    std::string coreShape = "PQ 28/20";
    std::string coreMaterial = "3C95";
    double ambientTemperature = 25;
    double coreLosses = 1.68;
    double expectedTemperature = 52;

    meanError += run_test_core_temperature(modelName, coreShape, coreMaterial, coreLosses, ambientTemperature,
                                           expectedTemperature);

    testCoreTemperatureAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core Temperature for " << coreMaterial << " with Model "
                  << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testCoreTemperatureAverageErrors[modelName].begin(),
                                 testCoreTemperatureAverageErrors[modelName].end()) /
                         testCoreTemperatureAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testCoreTemperatureMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

void test_core_temperature_miserable_43(CoreTemperatureModels modelName) {
    double meanError = 0;

    std::string coreShape = "PQ 26/20";
    std::string coreMaterial = "3C95";
    double ambientTemperature = 30;
    double coreLosses = 0.24;
    double expectedTemperature = 35;

    meanError += run_test_core_temperature(modelName, coreShape, coreMaterial, coreLosses, ambientTemperature,
                                           expectedTemperature);

    testCoreTemperatureAverageErrors[modelName].push_back(meanError);
    if (verboseTests) {
        std::cout << "Mean Error in Core Temperature for " << coreMaterial << " with Model "
                  << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testCoreTemperatureAverageErrors[modelName].begin(),
                                 testCoreTemperatureAverageErrors[modelName].end()) /
                         testCoreTemperatureAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        std::cout << "Current maximum for  " << magic_enum::enum_name(modelName) << ": "
                  << testCoreTemperatureMaximumErrors[modelName] * 100 << " %" << std::endl;
    }
}

SUITE(KazimierczukCoreTemperatureModel) {
    TEST(Test_Sotiris_47) {
        test_core_temperature_sotiris_47(CoreTemperatureModels::KAZIMIERCZUK);
    }
    TEST(Test_Sotiris_46) {
        test_core_temperature_sotiris_46(CoreTemperatureModels::KAZIMIERCZUK);
    }
    TEST(Test_Sotiris_40) {
        test_core_temperature_sotiris_40(CoreTemperatureModels::KAZIMIERCZUK);
    }
    TEST(Test_Sotiris_37) {
        test_core_temperature_sotiris_37(CoreTemperatureModels::KAZIMIERCZUK);
    }
    TEST(Test_Miserable_40) {
        test_core_temperature_miserable_40(CoreTemperatureModels::KAZIMIERCZUK);
    }
    TEST(Test_Miserable_43) {
        test_core_temperature_miserable_43(CoreTemperatureModels::KAZIMIERCZUK);
    }
}

SUITE(ManiktalaCoreTemperatureModel) {
    TEST(Test_Sotiris_47) {
        test_core_temperature_sotiris_47(CoreTemperatureModels::MANIKTALA);
    }
    TEST(Test_Sotiris_46) {
        test_core_temperature_sotiris_46(CoreTemperatureModels::MANIKTALA);
    }
    TEST(Test_Sotiris_40) {
        test_core_temperature_sotiris_40(CoreTemperatureModels::MANIKTALA);
    }
    TEST(Test_Sotiris_37) {
        test_core_temperature_sotiris_37(CoreTemperatureModels::MANIKTALA);
    }
    TEST(Test_Miserable_40) {
        test_core_temperature_miserable_40(CoreTemperatureModels::MANIKTALA);
    }
    TEST(Test_Miserable_43) {
        test_core_temperature_miserable_43(CoreTemperatureModels::MANIKTALA);
    }
}

SUITE(TdkCoreTemperatureModel) {
    TEST(Test_Sotiris_47) {
        test_core_temperature_sotiris_47(CoreTemperatureModels::TDK);
    }
    TEST(Test_Sotiris_46) {
        test_core_temperature_sotiris_46(CoreTemperatureModels::TDK);
    }
    TEST(Test_Sotiris_40) {
        test_core_temperature_sotiris_40(CoreTemperatureModels::TDK);
    }
    TEST(Test_Sotiris_37) {
        test_core_temperature_sotiris_37(CoreTemperatureModels::TDK);
    }
    TEST(Test_Miserable_40) {
        test_core_temperature_miserable_40(CoreTemperatureModels::TDK);
    }
    TEST(Test_Miserable_43) {
        test_core_temperature_miserable_43(CoreTemperatureModels::TDK);
    }
}

SUITE(DixonCoreTemperatureModel) {
    TEST(Test_Sotiris_47) {
        test_core_temperature_sotiris_47(CoreTemperatureModels::DIXON);
    }
    TEST(Test_Sotiris_46) {
        test_core_temperature_sotiris_46(CoreTemperatureModels::DIXON);
    }
    TEST(Test_Sotiris_40) {
        test_core_temperature_sotiris_40(CoreTemperatureModels::DIXON);
    }
    TEST(Test_Sotiris_37) {
        test_core_temperature_sotiris_37(CoreTemperatureModels::DIXON);
    }
    TEST(Test_Miserable_40) {
        test_core_temperature_miserable_40(CoreTemperatureModels::DIXON);
    }
    TEST(Test_Miserable_43) {
        test_core_temperature_miserable_43(CoreTemperatureModels::DIXON);
    }
}

SUITE(AmidonCoreTemperatureModel) {
    TEST(Test_Sotiris_47) {
        test_core_temperature_sotiris_47(CoreTemperatureModels::AMIDON);
    }
    TEST(Test_Sotiris_46) {
        test_core_temperature_sotiris_46(CoreTemperatureModels::AMIDON);
    }
    TEST(Test_Sotiris_40) {
        test_core_temperature_sotiris_40(CoreTemperatureModels::AMIDON);
    }
    TEST(Test_Sotiris_37) {
        test_core_temperature_sotiris_37(CoreTemperatureModels::AMIDON);
    }
    TEST(Test_Miserable_40) {
        test_core_temperature_miserable_40(CoreTemperatureModels::AMIDON);
    }
    TEST(Test_Miserable_43) {
        test_core_temperature_miserable_43(CoreTemperatureModels::AMIDON);
    }
}