#include "physical_models/CoreTemperature.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {  // Anonymous namespace to isolate test globals


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
    REQUIRE_THAT(calculatedTemperature, Catch::Matchers::WithinAbs(expectedCoreTemperature, expectedCoreTemperature * maximumAdmittedErrorTemperatureValue));

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

// ------------------------------------------------------------------
// Model-matrix tests: each fixture below used to be duplicated once per
// core-temperature model (5 near-identical 3-line TEST_CASEs per
// fixture, suffixed _2.._5). One TEST_CASE per fixture now GENERATEs
// over the models; INFO names the model so a failure is still
// attributable. The tag list keeps every per-model tag so tag-based
// selection still works.
// ------------------------------------------------------------------

#define ALL_CORE_TEMPERATURE_MODEL_TAGS                                                    \
    "[kazimierczuk-core-temperature-model][maniktala-core-temperature-model]"              \
    "[tdk-core-temperature-model][dixon-core-temperature-model]"                           \
    "[amidon-core-temperature-model]"

CoreTemperatureModels generate_all_core_temperature_models() {
    return GENERATE(CoreTemperatureModels::KAZIMIERCZUK,
                    CoreTemperatureModels::MANIKTALA,
                    CoreTemperatureModels::TDK,
                    CoreTemperatureModels::DIXON,
                    CoreTemperatureModels::AMIDON);
}

TEST_CASE("Test_Sotiris_47",
          "[physical-model][core-temperature]" ALL_CORE_TEMPERATURE_MODEL_TAGS "[smoke-test]") {
    auto modelName = generate_all_core_temperature_models();
    INFO("Model: " << magic_enum::enum_name(modelName));
    test_core_temperature_sotiris_47(modelName);
}

TEST_CASE("Test_Sotiris_46",
          "[physical-model][core-temperature]" ALL_CORE_TEMPERATURE_MODEL_TAGS "[smoke-test]") {
    auto modelName = generate_all_core_temperature_models();
    INFO("Model: " << magic_enum::enum_name(modelName));
    test_core_temperature_sotiris_46(modelName);
}

TEST_CASE("Test_Sotiris_40",
          "[physical-model][core-temperature]" ALL_CORE_TEMPERATURE_MODEL_TAGS "[smoke-test]") {
    auto modelName = generate_all_core_temperature_models();
    INFO("Model: " << magic_enum::enum_name(modelName));
    test_core_temperature_sotiris_40(modelName);
}

TEST_CASE("Test_Sotiris_37",
          "[physical-model][core-temperature]" ALL_CORE_TEMPERATURE_MODEL_TAGS "[smoke-test]") {
    auto modelName = generate_all_core_temperature_models();
    INFO("Model: " << magic_enum::enum_name(modelName));
    test_core_temperature_sotiris_37(modelName);
}

TEST_CASE("Test_Miserable_40",
          "[physical-model][core-temperature]" ALL_CORE_TEMPERATURE_MODEL_TAGS "[smoke-test]") {
    auto modelName = generate_all_core_temperature_models();
    INFO("Model: " << magic_enum::enum_name(modelName));
    test_core_temperature_miserable_40(modelName);
}

TEST_CASE("Test_Miserable_43",
          "[physical-model][core-temperature]" ALL_CORE_TEMPERATURE_MODEL_TAGS "[smoke-test]") {
    auto modelName = generate_all_core_temperature_models();
    INFO("Model: " << magic_enum::enum_name(modelName));
    test_core_temperature_miserable_43(modelName);
}

}  // namespace

// ABT #366/#362/#357: core-temperature models over the new families. These models read the
// core's effective volume / dimensions, which every new family now supplies through its own
// piece class — but nothing exercised that, and a thermal model quietly reading a zero or
// absent dimension would return a plausible-looking number. Smoke only: each family must
// produce a finite temperature above ambient that grows with dissipation. No pinned values —
// these are geometry classes MKF has no published thermal validation for.
TEST_CASE("Test_Core_Temperature_New_Families", "[physical-model][core-temperature][drum][drum-ring][drum-semishielded][molded]") {
    settings.reset();
    clear_databases();
    double ambientTemperature = 25;

    auto buildCustomCore = [](json shapeJson, const std::string& coreType, json coating = json()) {
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", coreType}, {"material", "3C90"}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        if (!coating.is_null()) {
            coreJson["functionalDescription"]["coating"] = coating;
        }
        Core core(coreJson);
        core.process_data();
        core.process_gap();
        return core;
    };

    json drumDimensions = {
        {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
        {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}}};
    json semishieldedDimensions = drumDimensions;
    semishieldedDimensions["J"] = {{"nominal", 0.0040}};
    semishieldedDimensions["K"] = {{"nominal", 0.0040}};
    semishieldedDimensions["L"] = {{"nominal", 0.0018}};

    std::vector<std::pair<std::string, Core>> cores;
    cores.emplace_back("drum", OpenMagneticsTesting::get_quick_core("DRH-14X20-4C", json::array(), 1, "3C90"));
    cores.emplace_back("drumRing", OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90"));
    cores.emplace_back("drumSemishielded", buildCustomCore(
        {{"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
         {"aliases", json::array()}, {"name", "LQS-like 4018"}, {"dimensions", semishieldedDimensions}},
        "pieceAndPlate",
        {{"type", "magneticEpoxy"}, {"thickness", 0.0001}, {"material", "Kool Mµ 26"}}));
    cores.emplace_back("molded", buildCustomCore(
        {{"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
         {"aliases", json::array()}, {"name", "MAPI-like 4020"},
         {"dimensions", {
             {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
             {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}},
        "closedShape"));

    auto coreTemperatureModel = CoreTemperatureModel::factory(CoreTemperatureModels::MANIKTALA);
    for (auto& [label, core] : cores) {
        double lowLossTemperature = coreTemperatureModel
            ->get_core_temperature(core, 0.01, ambientTemperature).get_maximum_temperature();
        double highLossTemperature = coreTemperatureModel
            ->get_core_temperature(core, 0.10, ambientTemperature).get_maximum_temperature();
        UNSCOPED_INFO(label << ": 10 mW -> " << lowLossTemperature << " C, 100 mW -> "
                      << highLossTemperature << " C");
        CHECK(std::isfinite(lowLossTemperature));
        CHECK(lowLossTemperature > ambientTemperature);
        CHECK(highLossTemperature > lowLossTemperature);
    }
    settings.reset();
}
