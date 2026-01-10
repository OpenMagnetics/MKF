#include "RandomUtils.h"
#include "support/Utils.h"
#include "physical_models/Reluctance.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    double maxError = 0.01;
    std::map<ReluctanceModels, double> maximumErrorReluctances = {
        {ReluctanceModels::ZHANG, 0.26},
        {ReluctanceModels::MUEHLETHALER, 0.42},
        {ReluctanceModels::EFFECTIVE_AREA, 0.42},
        {ReluctanceModels::EFFECTIVE_LENGTH, 0.42},
        {ReluctanceModels::PARTRIDGE, 0.32},
        {ReluctanceModels::STENGLEIN, 0.36},
        {ReluctanceModels::BALAKRISHNAN, 0.31},
        {ReluctanceModels::CLASSIC, 0.81},
    };
    std::map<ReluctanceModels, std::vector<double>> testAverageErrors = {
        {ReluctanceModels::ZHANG, {}},          {ReluctanceModels::MUEHLETHALER, {}},
        {ReluctanceModels::EFFECTIVE_AREA, {}}, {ReluctanceModels::EFFECTIVE_LENGTH, {}},
        {ReluctanceModels::PARTRIDGE, {}},      {ReluctanceModels::STENGLEIN, {}},
        {ReluctanceModels::BALAKRISHNAN, {}},   {ReluctanceModels::CLASSIC, {}},
    };
    std::map<ReluctanceModels, std::vector<std::pair<json, double>>> testErrorVersusC1 = {
        {ReluctanceModels::ZHANG, {}},          {ReluctanceModels::MUEHLETHALER, {}},
        {ReluctanceModels::EFFECTIVE_AREA, {}}, {ReluctanceModels::EFFECTIVE_LENGTH, {}},
        {ReluctanceModels::PARTRIDGE, {}},      {ReluctanceModels::STENGLEIN, {}},
        {ReluctanceModels::BALAKRISHNAN, {}},   {ReluctanceModels::CLASSIC, {}},
    };
    std::map<ReluctanceModels, double> maximumErrorEnergies = {
        {ReluctanceModels::ZHANG, 0.1},
        {ReluctanceModels::MUEHLETHALER, 0.23},
        {ReluctanceModels::EFFECTIVE_AREA, 0.13},
        {ReluctanceModels::EFFECTIVE_LENGTH, 0.13},
        {ReluctanceModels::PARTRIDGE, 0.1},
        {ReluctanceModels::STENGLEIN, 0.7},
        {ReluctanceModels::BALAKRISHNAN, 0.11},
        {ReluctanceModels::CLASSIC, 1},
    };

    double run_test_reluctance(ReluctanceModels modelName, std::string shapeName, json basicGapping,
                               double expectedReluctance) {
        double maximumErrorReluctance = maximumErrorReluctances[modelName];
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, basicGapping);
        auto gapping = core.get_functional_description().get_gapping();
        auto reluctanceModel = ReluctanceModel::factory(modelName);

        double calculatedReluctance = reluctanceModel->get_core_reluctance(core).get_core_reluctance();

        auto error = abs(expectedReluctance - calculatedReluctance) / expectedReluctance;
        std::pair<json, double> aux;
        json jsonAux;
        to_json(jsonAux, gapping[0]);
        aux.first = jsonAux;
        aux.second = error;
        testErrorVersusC1[modelName].push_back(aux);
        if (error > maximumErrorReluctance) {
            std::cout << "error " << error * 100 << " %" << std::endl;
        }
        REQUIRE_THAT(calculatedReluctance, Catch::Matchers::WithinAbs(expectedReluctance, expectedReluctance * maximumErrorReluctance));
        return error;
    }

    double run_test_energy(ReluctanceModels modelName, std::string shapeName, json basicGapping,
                           double expectedEnergy, int numberStacks = 1) {
        double maximumErrorEnergy = maximumErrorEnergies[modelName];
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, basicGapping, numberStacks);
        auto gapping = core.get_functional_description().get_gapping();
        double calculatedEnergy = 0;
        auto reluctanceModel = ReluctanceModel::factory(modelName);

        for (const auto& gap : gapping) {
            auto result = reluctanceModel->get_gap_reluctance(gap);
            calculatedEnergy += result.get_maximum_storable_magnetic_energy();

            REQUIRE(result.get_fringing_factor() >= 1);
            REQUIRE(result.get_maximum_storable_magnetic_energy() >= 0);
        }

        auto error = fabs(expectedEnergy - calculatedEnergy) / expectedEnergy;
        if (error > maximumErrorEnergy) {
            std::cout << "error " << error * 100 << " %" << std::endl;
        }
        REQUIRE_THAT(calculatedEnergy, Catch::Matchers::WithinAbs(expectedEnergy, expectedEnergy * maximumErrorEnergy));
        return error;
    }

    void test_reluctance_pq_28_20_grinded(ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {{{"gapLength", 0.4e-3}, {"expectedReluctance", 3446071}},
                                                            {{"gapLength", 0.5e-3}, {"expectedReluctance", 3233532}},
                                                            {{"gapLength", 0.7e-3}, {"expectedReluctance", 5514287}},
                                                            {{"gapLength", 1.08e-3}, {"expectedReluctance", 6871406}},
                                                            {{"gapLength", 1.65e-3}, {"expectedReluctance", 6982156}},
                                                            {{"gapLength", 0.305e-3}, {"expectedReluctance", 1736111}},
                                                            {{"gapLength", 0.305e-3}, {"expectedReluctance", 1736111}}};

        for (auto& test : tests) {
            meanError +=
                run_test_reluctance(modelName, "PQ 28/20", OpenMagneticsTesting::get_ground_gap(test["gapLength"]),
                                    test["expectedReluctance"]);
        }
        meanError /= tests.size();
        testAverageErrors[modelName].push_back(meanError);
        if (verboseTests) {
            std::cout << "Mean Error in Reluctance for Test_reluctance_PQ_28_20_Grinded with Model "
                      << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
            std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                  << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                         testAverageErrors[modelName].size() * 100
                  << " %" << std::endl;
        }
    }

    void test_reluctance_e42_21_20_spacer(ReluctanceModels modelName) {
        double meanError = 0;

        std::vector<std::map<std::string, double>> tests = {
            {{"gapLength", constants.residualGap}, {"expectedReluctance", 187891}},
            {{"gapLength", 0.0001}, {"expectedReluctance", 806451}},
            {{"gapLength", 0.00013}, {"expectedReluctance", 1035315}},
            {{"gapLength", 0.00015}, {"expectedReluctance", 1083841}},
            {{"gapLength", 0.00017}, {"expectedReluctance", 1358408}},
            {{"gapLength", 0.00020}, {"expectedReluctance", 1513877}},
            {{"gapLength", 0.0004}, {"expectedReluctance", 2441604}},
            {{"gapLength", 0.0005}, {"expectedReluctance", 3142238}},
            {{"gapLength", 0.001}, {"expectedReluctance", 4940440}},
        };

        for (auto& test : tests) {
            meanError +=
                run_test_reluctance(modelName, "E 42/21/20", OpenMagneticsTesting::get_spacer_gap(test["gapLength"]),
                                    test["expectedReluctance"]);
        }
        meanError /= tests.size();
        testAverageErrors[modelName].push_back(meanError);
        if (verboseTests) {
            std::cout << "Mean Error in Reluctance for Test_reluctance_E42_21_20_Spacer with Model "
                      << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
            std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                      << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                             testAverageErrors[modelName].size() * 100
                      << " %" << std::endl;
        }
    }

    void test_reluctance_etd_59_spacer(ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {
            {{"gapLength", 0.0001}, {"expectedReluctance", 565899}},
            {{"gapLength", 0.00013}, {"expectedReluctance", 698549}},
            {{"gapLength", 0.00015}, {"expectedReluctance", 752248}},
            {{"gapLength", 0.00017}, {"expectedReluctance", 905486}},
            {{"gapLength", 0.00020}, {"expectedReluctance", 1018686}},
            {{"gapLength", 0.0004}, {"expectedReluctance", 1610444}},
            {{"gapLength", 0.0005}, {"expectedReluctance", 2053962}},
            {{"gapLength", 0.001}, {"expectedReluctance", 3247502}},

        };

        for (auto& test : tests) {
            meanError +=
                run_test_reluctance(modelName, "ETD 59", OpenMagneticsTesting::get_spacer_gap(test["gapLength"]),
                                    test["expectedReluctance"]);
        }
        meanError /= tests.size();
        testAverageErrors[modelName].push_back(meanError);
        if (verboseTests) {
            std::cout << "Mean Error in Reluctance for Test_reluctance_ETD_59_Spacer with Model "
                      << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
            std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                      << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                             testAverageErrors[modelName].size() * 100
                      << " %" << std::endl;
        }
    }

    void test_reluctance_e_55_28_21_spacer(ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {
            {{"gapLength", 0.001}, {"expectedReluctance", 3091787}},
            {{"gapLength", 0.0015}, {"expectedReluctance", 4050632}},
            {{"gapLength", 0.002}, {"expectedReluctance", 5079365}},

        };

        for (auto& test : tests) {
            meanError +=
                run_test_reluctance(modelName, "E 55/28/21", OpenMagneticsTesting::get_spacer_gap(test["gapLength"]),
                                    test["expectedReluctance"]);
        }
        meanError /= tests.size();
        testAverageErrors[modelName].push_back(meanError);
        if (verboseTests) {
            std::cout << "Mean Error in Reluctance for Test_reluctance_E_55_28_21_Spacer with Model "
                      << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
            std::cout << "Current average for  " << magic_enum::enum_name(modelName) << ": "
                      << std::reduce(testAverageErrors[modelName].begin(), testAverageErrors[modelName].end()) /
                             testAverageErrors[modelName].size() * 100
                      << " %" << std::endl;

            // std::cout << "Current error versus C1 for  " << magic_enum::enum_name(modelName) << ": " << std::endl;
            // for (auto & ea: testErrorVersusC1[modelName]) {
            //     std::cout << ea.first << ": " << ea.second << std::endl;
            // }
        }
    }

    void test_energy_pq_40_40_grinded(ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {
            {{"gapLength", 0.002}, {"expectedEnergy", 0.02528}},

        };

        for (auto& test : tests) {
            meanError +=
                run_test_energy(modelName, "PQ 40/40", OpenMagneticsTesting::get_ground_gap(test["gapLength"]),
                                test["expectedEnergy"]);
        }
        meanError /= tests.size();
        if (verboseTests) {
            std::cout << "Mean Error in Energy for Test_reluctance_PQ_40_40_Grinded with Model "
                      << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
        }
    }

    TEST_CASE("Test_Zhang_PQ_28_20_Grinded", "[physical-model][reluctance][zhang-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::ZHANG);
    }

    TEST_CASE("Test_Zhang_E42_21_20_Spacer", "[physical-model][reluctance][zhang-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::ZHANG);
    }

    TEST_CASE("Test_Zhang_ETD_59_Spacer", "[physical-model][reluctance][zhang-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::ZHANG);
    }

    TEST_CASE("Test_Zhang_E_55_28_21_Spacer", "[physical-model][reluctance][zhang-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::ZHANG);
    }

    TEST_CASE("Test_Zhang_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][zhang-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.03255}};

        run_test_energy(ReluctanceModels::ZHANG, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Zhang_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][zhang-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.122}};

        run_test_energy(ReluctanceModels::ZHANG, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Zhang_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][zhang-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.2234}};

        run_test_energy(ReluctanceModels::ZHANG, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_Zhang_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][zhang-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.3245}};

        run_test_energy(ReluctanceModels::ZHANG, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_Zhang_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][zhang-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.092}};

        run_test_energy(ReluctanceModels::ZHANG, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_Muehlethaler_PQ_28_20_Grinded", "[physical-model][reluctance][muehlethaler-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::MUEHLETHALER);
    }

    TEST_CASE("Test_Muehlethaler_E42_21_20_Spacer", "[physical-model][reluctance][muehlethaler-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::MUEHLETHALER);
    }

    TEST_CASE("Test_Muehlethaler_ETD_59_Spacer", "[physical-model][reluctance][muehlethaler-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::MUEHLETHALER);
    }

    TEST_CASE("Test_Muehlethaler_E_55_28_21_Spacer", "[physical-model][reluctance][muehlethaler-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::MUEHLETHALER);
    }

    TEST_CASE("Test_Muehlethaler_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][muehlethaler-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.035}};

        run_test_energy(ReluctanceModels::MUEHLETHALER, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Muehlethaler_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][muehlethaler-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.142}};

        run_test_energy(ReluctanceModels::MUEHLETHALER, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Muehlethaler_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][muehlethaler-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.248}};

        run_test_energy(ReluctanceModels::MUEHLETHALER, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_Muehlethaler_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][muehlethaler-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.355}};

        run_test_energy(ReluctanceModels::MUEHLETHALER, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_Muehlethaler_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][muehlethaler-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.099}};

        run_test_energy(ReluctanceModels::MUEHLETHALER, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_EffectiveArea_PQ_28_20_Grinded", "[physical-model][reluctance][effective-area-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::EFFECTIVE_AREA);
    }

    TEST_CASE("Test_EffectiveArea_E42_21_20_Spacer", "[physical-model][reluctance][effective-area-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::EFFECTIVE_AREA);
    }

    TEST_CASE("Test_EffectiveArea_ETD_59_Spacer", "[physical-model][reluctance][effective-area-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::EFFECTIVE_AREA);
    }

    TEST_CASE("Test_EffectiveArea_E_55_28_21_Spacer", "[physical-model][reluctance][effective-area-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::EFFECTIVE_AREA);
    }

    TEST_CASE("Test_EffectiveArea_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][effective-area-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02872}};

        run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_EffectiveArea_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][effective-area-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1038}};

        run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_EffectiveArea_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][effective-area-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1945}};

        run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_EffectiveArea_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][effective-area-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.2852}};

        run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_EffectiveArea_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][effective-area-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.086}};

        run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_EffectiveLength_PQ_28_20_Grinded", "[physical-model][reluctance][effective-length-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::EFFECTIVE_LENGTH);
    }

    TEST_CASE("Test_EffectiveLength_E42_21_20_Spacer", "[physical-model][reluctance][effective-length-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::EFFECTIVE_LENGTH);
    }

    TEST_CASE("Test_EffectiveLength_ETD_59_Spacer", "[physical-model][reluctance][effective-length-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::EFFECTIVE_LENGTH);
    }

    TEST_CASE("Test_EffectiveLength_E_55_28_21_Spacer", "[physical-model][reluctance][effective-length-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::EFFECTIVE_LENGTH);
    }

    TEST_CASE("Test_EffectiveLength_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][effective-length-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02872}};

        run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_EffectiveLength_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][effective-length-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1038}};

        run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_EffectiveLength_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][effective-length-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1945}};

        run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_EffectiveLength_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][effective-length-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.2852}};

        run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_EffectiveLength_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][effective-length-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.086}};

        run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_Partridge_PQ_28_20_Grinded", "[physical-model][reluctance][partridge-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::PARTRIDGE);
    }

    TEST_CASE("Test_Partridge_E42_21_20_Spacer", "[physical-model][reluctance][partridge-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::PARTRIDGE);
    }

    TEST_CASE("Test_Partridge_ETD_59_Spacer", "[physical-model][reluctance][partridge-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::PARTRIDGE);
    }

    TEST_CASE("Test_Partridge_E_55_28_21_Spacer", "[physical-model][reluctance][partridge-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::PARTRIDGE);
    }

    TEST_CASE("Test_Partridge_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][partridge-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.033}};

        run_test_energy(ReluctanceModels::PARTRIDGE, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Partridge_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][partridge-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.12}};

        run_test_energy(ReluctanceModels::PARTRIDGE, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Partridge_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][partridge-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.216}};

        run_test_energy(ReluctanceModels::PARTRIDGE, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_Partridge_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][partridge-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.308}};

        run_test_energy(ReluctanceModels::PARTRIDGE, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_Partridge_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][partridge-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.095}};

        run_test_energy(ReluctanceModels::PARTRIDGE, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_Stenglein_PQ_28_20_Grinded", "[physical-model][reluctance][stenglein-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::STENGLEIN);
    }

    TEST_CASE("Test_Stenglein_E42_21_20_Spacer", "[physical-model][reluctance][stenglein-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::STENGLEIN);
    }

    TEST_CASE("Test_Stenglein_ETD_59_Spacer", "[physical-model][reluctance][stenglein-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::STENGLEIN);
    }

    TEST_CASE("Test_Stenglein_E_55_28_21_Spacer", "[physical-model][reluctance][stenglein-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::STENGLEIN);
    }

    TEST_CASE("Test_Stenglein_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][stenglein-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02528}};

        run_test_energy(ReluctanceModels::STENGLEIN, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Stenglein_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][stenglein-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.09}};

        run_test_energy(ReluctanceModels::STENGLEIN, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Stenglein_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][stenglein-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1629}};

        run_test_energy(ReluctanceModels::STENGLEIN, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_Stenglein_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][stenglein-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.24}};

        run_test_energy(ReluctanceModels::STENGLEIN, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_Stenglein_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][stenglein-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.165}};

        run_test_energy(ReluctanceModels::STENGLEIN, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_Balakrishnan_PQ_28_20_Grinded", "[physical-model][reluctance][balakrishnan-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::BALAKRISHNAN);
    }

    TEST_CASE("Test_Balakrishnan_E42_21_20_Spacer", "[physical-model][reluctance][balakrishnan-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::BALAKRISHNAN);
    }

    TEST_CASE("Test_Balakrishnan_ETD_59_Spacer", "[physical-model][reluctance][balakrishnan-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::BALAKRISHNAN);
    }

    TEST_CASE("Test_Balakrishnan_E_55_28_21_Spacer", "[physical-model][reluctance][balakrishnan-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::BALAKRISHNAN);
    }

    TEST_CASE("Test_Balakrishnan_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][balakrishnan-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.030}};

        run_test_energy(ReluctanceModels::BALAKRISHNAN, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Balakrishnan_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][balakrishnan-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.106}};

        run_test_energy(ReluctanceModels::BALAKRISHNAN, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Balakrishnan_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][balakrishnan-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.212}};

        run_test_energy(ReluctanceModels::BALAKRISHNAN, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_Balakrishnan_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][balakrishnan-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.319}};

        run_test_energy(ReluctanceModels::BALAKRISHNAN, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_Balakrishnan_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][balakrishnan-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.087}};

        run_test_energy(ReluctanceModels::BALAKRISHNAN, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_Classic_PQ_28_20_Grinded", "[physical-model][reluctance][classic-reluctance-model]") {
        test_reluctance_pq_28_20_grinded(ReluctanceModels::CLASSIC);
    }

    TEST_CASE("Test_Classic_E42_21_20_Spacer", "[physical-model][reluctance][classic-reluctance-model]") {
        test_reluctance_e42_21_20_spacer(ReluctanceModels::CLASSIC);
    }

    TEST_CASE("Test_Classic_ETD_59_Spacer", "[physical-model][reluctance][classic-reluctance-model]") {
        test_reluctance_etd_59_spacer(ReluctanceModels::CLASSIC);
    }

    TEST_CASE("Test_Classic_E_55_28_21_Spacer", "[physical-model][reluctance][classic-reluctance-model]") {
        test_reluctance_e_55_28_21_spacer(ReluctanceModels::CLASSIC);
    }

    TEST_CASE("Test_Classic_Energy_PQ_40_40_Grinded", "[physical-model][reluctance][magnetic-energy][classic-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02528}};

        run_test_energy(ReluctanceModels::CLASSIC, "PQ 40/40",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Classic_Energy_E_80_38_20_Grinded", "[physical-model][reluctance][magnetic-energy][classic-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.09}};

        run_test_energy(ReluctanceModels::CLASSIC, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
    }

    TEST_CASE("Test_Classic_Energy_E_80_38_20_2_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][classic-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1629}};

        run_test_energy(ReluctanceModels::CLASSIC, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
    }

    TEST_CASE("Test_Classic_Energy_E_80_38_20_3_Stacks_Grinded", "[physical-model][reluctance][magnetic-energy][classic-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.24}};

        run_test_energy(ReluctanceModels::CLASSIC, "E 80/38/20",
                        OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
    }

    TEST_CASE("Test_Classic_Energy_PQ_40_40_Distributed", "[physical-model][reluctance][magnetic-energy][classic-reluctance-model]") {
        std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.076}};

        run_test_energy(ReluctanceModels::CLASSIC, "PQ 40/40",
                        OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
    }

    TEST_CASE("Test_Gap_By_Fringing_Factor", "[physical-model][reluctance][gap]") {
        for (size_t i = 0; i < 100; ++i)
        {
            double randomPercentage = double(OpenMagnetics::TestUtils::randomInt64(1, 50 + 1 - 1)) * 1e-2;
            auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", OpenMagneticsTesting::get_residual_gap());
            auto centralColumns = core.find_columns_by_type(ColumnType::CENTRAL);
            double expectedGap = centralColumns[0].get_height() * randomPercentage;
            core = OpenMagneticsTesting::get_quick_core("E 42/21/20", OpenMagneticsTesting::get_ground_gap(expectedGap));
            auto reluctanceModel = ReluctanceModel::factory(ReluctanceModels::ZHANG);
            auto fringingFactor = reluctanceModel->get_core_reluctance(core).get_maximum_fringing_factor().value();

            double gap = reluctanceModel->get_gapping_by_fringing_factor(core, fringingFactor);
            REQUIRE_THAT(expectedGap, Catch::Matchers::WithinAbs(gap, expectedGap * maxError));
        }
    }

    TEST_CASE("Test_Reluctance_3C96", "[physical-model][reluctance]") {
        double dcCurrent = 10;
        double ambientTemperature = 50;
        double dutyCycle = 0.75;
        double peakToPeak = 10;
        double frequency = 100000;
        double magnetizingInductance = 10e-6;

        OpenMagnetics::Inputs inputs = OpenMagnetics::Inputs::create_quick_operating_point(
            frequency, magnetizingInductance, ambientTemperature, WaveformLabel::RECTANGULAR, peakToPeak,
            dutyCycle, dcCurrent);
        auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", OpenMagneticsTesting::get_residual_gap());
        auto reluctanceModel = ReluctanceModel::factory(ReluctanceModels::ZHANG);

        double calculatedReluctanceAt50 =
            reluctanceModel->get_ungapped_core_reluctance(core, inputs.get_mutable_operating_points()[0]);

        ambientTemperature = 100;
        inputs = OpenMagnetics::Inputs::create_quick_operating_point(
            frequency, magnetizingInductance, ambientTemperature, WaveformLabel::RECTANGULAR, peakToPeak,
            dutyCycle, dcCurrent);
        double calculatedReluctanceAt100 =
            reluctanceModel->get_ungapped_core_reluctance(core, inputs.get_mutable_operating_points()[0]);

        ambientTemperature = 150;
        inputs = OpenMagnetics::Inputs::create_quick_operating_point(
            frequency, magnetizingInductance, ambientTemperature, WaveformLabel::RECTANGULAR, peakToPeak,
            dutyCycle, dcCurrent);
        double calculatedReluctanceAt150 =
            reluctanceModel->get_ungapped_core_reluctance(core, inputs.get_mutable_operating_points()[0]);

        ambientTemperature = 200;
        inputs = OpenMagnetics::Inputs::create_quick_operating_point(
            frequency, magnetizingInductance, ambientTemperature, WaveformLabel::RECTANGULAR, peakToPeak,
            dutyCycle, dcCurrent);
        double calculatedReluctanceAt200 =
            reluctanceModel->get_ungapped_core_reluctance(core, inputs.get_mutable_operating_points()[0]);

        REQUIRE(calculatedReluctanceAt50 > calculatedReluctanceAt100);
        REQUIRE(calculatedReluctanceAt150 > calculatedReluctanceAt100);
        REQUIRE(calculatedReluctanceAt150 > calculatedReluctanceAt200);
    }

    TEST_CASE("Test_Web_0", "[physical-model][reluctance][gap][bug]") {
        std::string coreGapData = R"({"area":0.000123,"coordinates":[0,0.0005,0],"distanceClosestNormalSurface":0.014098,"distanceClosestParallelSurface":0.0088,"length":0.000005,"sectionDimensions":[0.0125,0.0125],"shape":"round","type":"subtractive"})";
        std::string modelNameString = "Zhang";
        std::string modelNameStringUpper = modelNameString;
        std::transform(modelNameStringUpper.begin(), modelNameStringUpper.end(), modelNameStringUpper.begin(), ::toupper);
        auto modelName = magic_enum::enum_cast<ReluctanceModels>(modelNameStringUpper);

        auto reluctanceModel = ReluctanceModel::factory(modelName.value());

        CoreGap coreGap(json::parse(coreGapData));

        auto coreGapResult = reluctanceModel->get_gap_reluctance(coreGap);
    }

}  // namespace
