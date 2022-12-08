#include <UnitTest++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include "Reluctance.h"
#include "Constants.h"


SUITE(Reluctance)
{
    std::map<OpenMagnetics::ReluctanceModels, double> maximumErrors = {
        {OpenMagnetics::ReluctanceModels::ZHANG, 0.26},
        {OpenMagnetics::ReluctanceModels::MUEHLETHALER, 0.41},
        {OpenMagnetics::ReluctanceModels::EFFECTIVE_AREA, 0.42},
        {OpenMagnetics::ReluctanceModels::MCLYMAN, 0.32},
    };

    OpenMagnetics::Core get_core(std::string shapeName, json basicGapping) {
        auto coreJson = json();

        coreJson["functionalDescription"] = json();
        coreJson["functionalDescription"]["name"] = "GapReluctanceTest";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N87";
        coreJson["functionalDescription"]["shape"] = shapeName;
        coreJson["functionalDescription"]["gapping"] = basicGapping;
        coreJson["functionalDescription"]["numberStacks"] = 1;

        OpenMagnetics::Core core(coreJson);

        return core;
    }

    json get_grinded_gap(double gapLength) {
        auto constants = OpenMagnetics::Constants();
        auto basicGapping = json::array();
        auto basicCentralGap = json();
        basicCentralGap["type"] = "subtractive";
        basicCentralGap["length"] = gapLength;
        auto basicLateralGap = json();
        basicLateralGap["type"] = "residual";
        basicLateralGap["length"] = constants.residualGap;
        basicGapping.push_back(basicCentralGap);
        basicGapping.push_back(basicLateralGap);
        basicGapping.push_back(basicLateralGap);

        return basicGapping;
    }

    json get_spacer_gap(double gapLength) {
        auto basicGapping = json::array();
        auto basicSpacerGap = json();
        basicSpacerGap["type"] = "additive";
        basicSpacerGap["length"] = gapLength;
        basicGapping.push_back(basicSpacerGap);
        basicGapping.push_back(basicSpacerGap);
        basicGapping.push_back(basicSpacerGap);

        return basicGapping;
    }

    double run_test(OpenMagnetics::ReluctanceModels modelName, std::string shapeName, json basicGapping, double expectedReluctance) {

        double maximumError = maximumErrors[modelName];
        auto core = get_core(shapeName, basicGapping);
        auto gapping = core.get_functional_description().get_gapping();
        double calculatedReluctance = 0;
        double calculatedCentralReluctance = 0;
        double calculatedLateralReluctance = 0;
        auto reluctanceZhangModel = OpenMagnetics::ReluctanceModel::factory(modelName);
        auto coreReluctance = reluctanceZhangModel->get_core_reluctance(core);

        for (const auto& gap: gapping) {
            auto result = reluctanceZhangModel->get_gap_reluctance(gap);
            auto gapColumn = core.find_closest_column_by_coordinates(*gap.get_coordinates());
            if (gapColumn.get_type() == OpenMagnetics::ColumnType::LATERAL) {
                calculatedLateralReluctance += 1 / result["reluctance"];
            }
            else {
                calculatedCentralReluctance += result["reluctance"];
            }
            if (result["fringing_factor"] < 1) {
                std::cout << "fringing_factor " << result["fringing_factor"] << std::endl;
            }
            CHECK(result["fringing_factor"] >= 1);

        }
        calculatedReluctance = coreReluctance + calculatedCentralReluctance + 1 / calculatedLateralReluctance;

        auto error = abs(expectedReluctance - calculatedReluctance) / expectedReluctance;
        if (error > maximumError) {
            std::cout << "error " << error * 100 << " %" << std::endl;

        }
        CHECK_CLOSE(calculatedReluctance, expectedReluctance, expectedReluctance * maximumError);
        return error;
    }

    void test_pq_28_20_grinded(OpenMagnetics::ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {
        {
            {"gapLength", 0.4e-3},
            {"expectedReluctance", 3446071}
        },
        {
            {"gapLength", 0.5e-3},
            {"expectedReluctance", 3233532}
        },
        {
            {"gapLength", 0.7e-3},
            {"expectedReluctance", 5514287}
        },
        {
            {"gapLength", 1.08e-3},
            {"expectedReluctance", 6871406}
        },
        {
            {"gapLength", 1.65e-3},
            {"expectedReluctance", 6982156}
        },
        {
            {"gapLength", 0.305e-3},
            {"expectedReluctance", 1736111}
        },
        {
            {"gapLength", 0.305e-3},
            {"expectedReluctance", 1736111}
        }};

        for (auto& test : tests) {
            meanError += run_test(modelName, "PQ 28/20", get_grinded_gap(test["gapLength"]), test["expectedReluctance"]);
        }
        meanError /= tests.size();
        std::cout << "Mean Error for Test_PQ_28_20_Grinded with Model " << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
    }

    void test_e42_21_20_spacer(OpenMagnetics::ReluctanceModels modelName) {
        auto constants = OpenMagnetics::Constants();
        double meanError = 0;

        std::vector<std::map<std::string, double>> tests = {
            {
                {"gapLength", constants.residualGap},
                {"expectedReluctance", 187891}
            },
            {
                {"gapLength", 0.0001},
                {"expectedReluctance", 806451}
            },
            {
                {"gapLength", 0.00013},
                {"expectedReluctance", 1035315}
            },
            {
                {"gapLength", 0.00015},
                {"expectedReluctance", 1083841}
            },
            {
                {"gapLength", 0.00017},
                {"expectedReluctance", 1358408}
            },
            {
                {"gapLength", 0.00020},
                {"expectedReluctance", 1513877}
            },
            {
                {"gapLength", 0.0004},
                {"expectedReluctance", 2441604}
            },
            {
                {"gapLength", 0.0005},
                {"expectedReluctance", 3142238}
            },
            {
                {"gapLength", 0.001},
                {"expectedReluctance", 4940440}
            },
        };

        for (auto& test : tests) {
            meanError += run_test(modelName, "E 42/21/20", get_spacer_gap(test["gapLength"]), test["expectedReluctance"]);
        }
        meanError /= tests.size();
        std::cout << "Mean Error for Test_E42_21_20_Spacer with Model " << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
    
    }

    void test_etd_59_spacer(OpenMagnetics::ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {
            {
                {"gapLength", 0.0001},
                {"expectedReluctance", 565899}
            },
            {
                {"gapLength", 0.00013},
                {"expectedReluctance", 698549}
            },
            {
                {"gapLength", 0.00015},
                {"expectedReluctance", 752248}
            },
            {
                {"gapLength", 0.00017},
                {"expectedReluctance", 905486}
            },
            {
                {"gapLength", 0.00020},
                {"expectedReluctance", 1018686}
            },
            {
                {"gapLength", 0.0004},
                {"expectedReluctance", 1610444}
            },
            {
                {"gapLength", 0.0005},
                {"expectedReluctance", 2053962}
            },
            {
                {"gapLength", 0.001},
                {"expectedReluctance", 3247502}
            },

        };

        for (auto& test : tests) {
            meanError += run_test(modelName, "ETD 59", get_spacer_gap(test["gapLength"]), test["expectedReluctance"]);
        }
        meanError /= tests.size();
        std::cout << "Mean Error for Test_ETD_59_Spacer with Model " << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;    
    }

    void test_e_55_28_21_spacer(OpenMagnetics::ReluctanceModels modelName) {
        double meanError = 0;
        std::vector<std::map<std::string, double>> tests = {
            {
                {"gapLength", 0.001},
                {"expectedReluctance", 3091787}
            },
            {
                {"gapLength", 0.0015},
                {"expectedReluctance", 4050632}
            },
            {
                {"gapLength", 0.002},
                {"expectedReluctance", 5079365}
            },

        };

        for (auto& test : tests) {
            meanError += run_test(modelName, "E 55/28/21", get_spacer_gap(test["gapLength"]), test["expectedReluctance"]);
        }
        meanError /= tests.size();
        std::cout << "Mean Error for Test_E_55_28_21_Spacer with Model " << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
    
    }

    SUITE(ZhangModel)
    {
        TEST(Test_PQ_28_20_Grinded)
        {
            test_pq_28_20_grinded(OpenMagnetics::ReluctanceModels::ZHANG);
        }

        TEST(Test_E42_21_20_Spacer)
        {
            test_e42_21_20_spacer(OpenMagnetics::ReluctanceModels::ZHANG);
        }

        TEST(Test_ETD_59_Spacer)
        {
            test_etd_59_spacer(OpenMagnetics::ReluctanceModels::ZHANG);
        }

        TEST(Test_E_55_28_21_Spacer)
        {
            test_e_55_28_21_spacer(OpenMagnetics::ReluctanceModels::ZHANG);
        }
    }

    SUITE(MuehlethalerModel)
    {
        TEST(Test_PQ_28_20_Grinded)
        {
            test_pq_28_20_grinded(OpenMagnetics::ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_E42_21_20_Spacer)
        {
            test_e42_21_20_spacer(OpenMagnetics::ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_ETD_59_Spacer)
        {
            test_etd_59_spacer(OpenMagnetics::ReluctanceModels::MUEHLETHALER);
        }

        TEST(Test_E_55_28_21_Spacer)
        {
            test_e_55_28_21_spacer(OpenMagnetics::ReluctanceModels::MUEHLETHALER);
        }
    }

    SUITE(EffectiveAreaModel)
    {
        TEST(Test_PQ_28_20_Grinded)
        {
            test_pq_28_20_grinded(OpenMagnetics::ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_E42_21_20_Spacer)
        {
            test_e42_21_20_spacer(OpenMagnetics::ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_ETD_59_Spacer)
        {
            test_etd_59_spacer(OpenMagnetics::ReluctanceModels::EFFECTIVE_AREA);
        }

        TEST(Test_E_55_28_21_Spacer)
        {
            test_e_55_28_21_spacer(OpenMagnetics::ReluctanceModels::EFFECTIVE_AREA);
        }
    }

    SUITE(McLymanModel)
    {
        TEST(Test_PQ_28_20_Grinded)
        {
            test_pq_28_20_grinded(OpenMagnetics::ReluctanceModels::MCLYMAN);
        }

        TEST(Test_E42_21_20_Spacer)
        {
            test_e42_21_20_spacer(OpenMagnetics::ReluctanceModels::MCLYMAN);
        }

        TEST(Test_ETD_59_Spacer)
        {
            test_etd_59_spacer(OpenMagnetics::ReluctanceModels::MCLYMAN);
        }

        TEST(Test_E_55_28_21_Spacer)
        {
            test_e_55_28_21_spacer(OpenMagnetics::ReluctanceModels::MCLYMAN);
        }
    }
}
