#include "support/Utils.h"
#include "physical_models/Reluctance.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(Reluctance) {
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
        CHECK_CLOSE(calculatedReluctance, expectedReluctance, expectedReluctance * maximumErrorReluctance);
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

            CHECK(result.get_fringing_factor() >= 1);
            CHECK(result.get_maximum_storable_magnetic_energy() >= 0);
        }

        auto error = fabs(expectedEnergy - calculatedEnergy) / expectedEnergy;
        if (error > maximumErrorEnergy) {
            std::cout << "error " << error * 100 << " %" << std::endl;
        }
        CHECK_CLOSE(calculatedEnergy, expectedEnergy, expectedEnergy * maximumErrorEnergy);
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

    SUITE(ZhangModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::ZHANG);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::ZHANG);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::ZHANG);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::ZHANG);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.03255}};

            run_test_energy(ReluctanceModels::ZHANG, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.122}};

            run_test_energy(ReluctanceModels::ZHANG, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.2234}};

            run_test_energy(ReluctanceModels::ZHANG, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.3245}};

            run_test_energy(ReluctanceModels::ZHANG, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.092}};

            run_test_energy(ReluctanceModels::ZHANG, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(MuehlethalerModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.035}};

            run_test_energy(ReluctanceModels::MUEHLETHALER, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.142}};

            run_test_energy(ReluctanceModels::MUEHLETHALER, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.248}};

            run_test_energy(ReluctanceModels::MUEHLETHALER, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.355}};

            run_test_energy(ReluctanceModels::MUEHLETHALER, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.099}};

            run_test_energy(ReluctanceModels::MUEHLETHALER, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(EffectiveAreaModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02872}};

            run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1038}};

            run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1945}};

            run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.2852}};

            run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.086}};

            run_test_energy(ReluctanceModels::EFFECTIVE_AREA, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(EffectiveLengthModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::EFFECTIVE_LENGTH);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::EFFECTIVE_LENGTH);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::EFFECTIVE_LENGTH);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::EFFECTIVE_LENGTH);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02872}};

            run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1038}};

            run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1945}};

            run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.2852}};

            run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.086}};

            run_test_energy(ReluctanceModels::EFFECTIVE_LENGTH, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(PartridgeModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::PARTRIDGE);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::PARTRIDGE);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::PARTRIDGE);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::PARTRIDGE);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.033}};

            run_test_energy(ReluctanceModels::PARTRIDGE, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.12}};

            run_test_energy(ReluctanceModels::PARTRIDGE, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.216}};

            run_test_energy(ReluctanceModels::PARTRIDGE, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.308}};

            run_test_energy(ReluctanceModels::PARTRIDGE, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.095}};

            run_test_energy(ReluctanceModels::PARTRIDGE, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(StengleinModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::STENGLEIN);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::STENGLEIN);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::STENGLEIN);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::STENGLEIN);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02528}};

            run_test_energy(ReluctanceModels::STENGLEIN, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.09}};

            run_test_energy(ReluctanceModels::STENGLEIN, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1629}};

            run_test_energy(ReluctanceModels::STENGLEIN, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.24}};

            run_test_energy(ReluctanceModels::STENGLEIN, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.165}};

            run_test_energy(ReluctanceModels::STENGLEIN, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(BalakrishnanModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::BALAKRISHNAN);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::BALAKRISHNAN);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::BALAKRISHNAN);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::BALAKRISHNAN);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.030}};

            run_test_energy(ReluctanceModels::BALAKRISHNAN, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.106}};

            run_test_energy(ReluctanceModels::BALAKRISHNAN, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.212}};

            run_test_energy(ReluctanceModels::BALAKRISHNAN, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.319}};

            run_test_energy(ReluctanceModels::BALAKRISHNAN, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.087}};

            run_test_energy(ReluctanceModels::BALAKRISHNAN, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    SUITE(ClassicModel) {
        TEST(Test_PQ_28_20_Grinded) {
            test_reluctance_pq_28_20_grinded(ReluctanceModels::CLASSIC);
        }

        TEST(Test_E42_21_20_Spacer) {
            test_reluctance_e42_21_20_spacer(ReluctanceModels::CLASSIC);
        }

        TEST(Test_ETD_59_Spacer) {
            test_reluctance_etd_59_spacer(ReluctanceModels::CLASSIC);
        }

        TEST(Test_E_55_28_21_Spacer) {
            test_reluctance_e_55_28_21_spacer(ReluctanceModels::CLASSIC);
        }

        TEST(Test_Energy_PQ_40_40_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.02528}};

            run_test_energy(ReluctanceModels::CLASSIC, "PQ 40/40",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.09}};

            run_test_energy(ReluctanceModels::CLASSIC, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"]);
        }

        TEST(Test_Energy_E_80_38_20_2_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.1629}};

            run_test_energy(ReluctanceModels::CLASSIC, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 2);
        }

        TEST(Test_Energy_E_80_38_20_3_Stacks_Grinded) {
            std::map<std::string, double> test = {{"gapLength", 0.003}, {"expectedEnergy", 0.24}};

            run_test_energy(ReluctanceModels::CLASSIC, "E 80/38/20",
                            OpenMagneticsTesting::get_ground_gap(test["gapLength"]), test["expectedEnergy"], 3);
        }

        TEST(Test_Energy_PQ_40_40_Distributed) {
            std::map<std::string, double> test = {{"gapLength", 0.002}, {"expectedEnergy", 0.076}};

            run_test_energy(ReluctanceModels::CLASSIC, "PQ 40/40",
                            OpenMagneticsTesting::get_distributed_gap(test["gapLength"], 3), test["expectedEnergy"], 1);
        }
    }

    TEST(Test_Gap_By_Fringing_Factor) {
        srand (time(NULL));
        for (size_t i = 0; i < 100; ++i)
        {
            double randomPercentage = double(std::rand() % 50 + 1L) * 1e-2;
            auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", OpenMagneticsTesting::get_residual_gap());
            auto centralColumns = core.find_columns_by_type(ColumnType::CENTRAL);
            double expectedGap = centralColumns[0].get_height() * randomPercentage;
            core = OpenMagneticsTesting::get_quick_core("E 42/21/20", OpenMagneticsTesting::get_ground_gap(expectedGap));
            auto reluctanceModel = ReluctanceModel::factory(ReluctanceModels::ZHANG);
            auto fringingFactor = reluctanceModel->get_core_reluctance(core).get_maximum_fringing_factor().value();

            double gap = reluctanceModel->get_gapping_by_fringing_factor(core, fringingFactor);
            CHECK_CLOSE(expectedGap, gap, expectedGap * maxError);
        }
    }
}

SUITE(ReluctanceUngappedCore) {
    TEST(Test_Reluctance_3C96) {
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

        CHECK(calculatedReluctanceAt50 > calculatedReluctanceAt100);
        CHECK(calculatedReluctanceAt150 > calculatedReluctanceAt100);
        CHECK(calculatedReluctanceAt150 > calculatedReluctanceAt200);
    }
}


SUITE(WebReluctance) {
    TEST(Test_Web_0) {
        std::string coreGapData = R"({"area":0.000123,"coordinates":[0,0.0005,0],"distanceClosestNormalSurface":0.014098,"distanceClosestParallelSurface":0.0088,"length":0.000005,"sectionDimensions":[0.0125,0.0125],"shape":"round","type":"subtractive"})";
        std::string modelNameString = "Zhang";
        std::string modelNameStringUpper = modelNameString;
        std::transform(modelNameStringUpper.begin(), modelNameStringUpper.end(), modelNameStringUpper.begin(), ::toupper);
        auto modelName = magic_enum::enum_cast<ReluctanceModels>(modelNameStringUpper);

        auto reluctanceModel = ReluctanceModel::factory(modelName.value());

        CoreGap coreGap(json::parse(coreGapData));

        auto coreGapResult = reluctanceModel->get_gap_reluctance(coreGap);
    }
}
