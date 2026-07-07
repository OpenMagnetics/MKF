#include "RandomUtils.h"
#include "support/Utils.h"
#include "physical_models/Reluctance.h"
#include "TestingUtils.h"

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
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
        {ReluctanceModels::ZHANG, 0.27},
        {ReluctanceModels::MUEHLETHALER, 0.42},
        {ReluctanceModels::EFFECTIVE_AREA, 0.42},
        {ReluctanceModels::EFFECTIVE_LENGTH, 0.42},
        {ReluctanceModels::PARTRIDGE, 0.32},
        {ReluctanceModels::STENGLEIN, 0.36},
        {ReluctanceModels::BALAKRISHNAN, 0.32},
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

    // NOTE (July 2026 re-baseline): the expectedEnergy constants passed to run_test_energy were
    // reduced (by a per-model factor F^2) after ReluctanceModel::get_gap_maximum_storable_energy
    // was fixed to divide by the fringing factor F instead of multiplying by it. Gap stored energy
    // at core saturation is E = 0.5*Phi^2*R_gap with R_gap = lg/(mu0*A*F), i.e. E is proportional to
    // 1/F — the old *F over-counted by F^2. This brought Reluctance.h in line with the already-fixed
    // MagneticEnergy::get_gap_maximum_magnetic_energy. Values are the measured /F outputs.
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

    // Commented out: unused function causing -Wunused-function warning
    // void test_energy_pq_40_40_grinded(ReluctanceModels modelName) {
    //     double meanError = 0;
    //     std::vector<std::map<std::string, double>> tests = {
    //         {{"gapLength", 0.002}, {"expectedEnergy", 0.02528}},
    //
    //     };
    //
    //     for (auto& test : tests) {
    //         meanError +=
    //             run_test_energy(modelName, "PQ 40/40", OpenMagneticsTesting::get_ground_gap(test["gapLength"]),
    //                             test["expectedEnergy"]);
    //     }
    //     meanError /= tests.size();
    //     if (verboseTests) {
    //         std::cout << "Mean Error in Energy for Test_reluctance_PQ_40_40_Grinded with Model "
    //                   << magic_enum::enum_name(modelName) << ": " << meanError * 100 << " %" << std::endl;
    //     }
    // }

    // ------------------------------------------------------------------
    // Model-matrix tests: each fixture below used to be duplicated once
    // per reluctance model (8 near-identical 3-line TEST_CASEs per
    // fixture). One TEST_CASE per fixture now GENERATEs over the models;
    // INFO names the model so a failure is still attributable. The tag
    // list keeps every per-model tag so tag-based selection still works.
    // ------------------------------------------------------------------

#define ALL_RELUCTANCE_MODEL_TAGS                                                          \
    "[zhang-reluctance-model][muehlethaler-reluctance-model]"                              \
    "[effective-area-reluctance-model][effective-length-reluctance-model]"                 \
    "[partridge-reluctance-model][stenglein-reluctance-model]"                             \
    "[balakrishnan-reluctance-model][classic-reluctance-model]"

    ReluctanceModels generate_all_reluctance_models() {
        return GENERATE(ReluctanceModels::ZHANG,
                        ReluctanceModels::MUEHLETHALER,
                        ReluctanceModels::EFFECTIVE_AREA,
                        ReluctanceModels::EFFECTIVE_LENGTH,
                        ReluctanceModels::PARTRIDGE,
                        ReluctanceModels::STENGLEIN,
                        ReluctanceModels::BALAKRISHNAN,
                        ReluctanceModels::CLASSIC);
    }

    TEST_CASE("Test_Reluctance_PQ_28_20_Grinded",
              "[physical-model][reluctance]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto modelName = generate_all_reluctance_models();
        INFO("Model: " << magic_enum::enum_name(modelName));
        test_reluctance_pq_28_20_grinded(modelName);
    }

    TEST_CASE("Test_Reluctance_E42_21_20_Spacer",
              "[physical-model][reluctance]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto modelName = generate_all_reluctance_models();
        INFO("Model: " << magic_enum::enum_name(modelName));
        test_reluctance_e42_21_20_spacer(modelName);
    }

    TEST_CASE("Test_Reluctance_ETD_59_Spacer",
              "[physical-model][reluctance]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto modelName = generate_all_reluctance_models();
        INFO("Model: " << magic_enum::enum_name(modelName));
        test_reluctance_etd_59_spacer(modelName);
    }

    TEST_CASE("Test_Reluctance_E_55_28_21_Spacer",
              "[physical-model][reluctance]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto modelName = generate_all_reluctance_models();
        INFO("Model: " << magic_enum::enum_name(modelName));
        test_reluctance_e_55_28_21_spacer(modelName);
    }

    // Energy fixtures: expected energy depends on (model, fixture), so the
    // generator is a (model, expectedEnergy) table. Values are unchanged
    // from the per-model TEST_CASEs they replace.

    TEST_CASE("Test_Reluctance_Energy_PQ_40_40_Grinded",
              "[physical-model][reluctance][magnetic-energy]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto [modelName, expectedEnergy] = GENERATE(table<ReluctanceModels, double>({
            {ReluctanceModels::ZHANG, 0.015136975},
            {ReluctanceModels::MUEHLETHALER, 0.014135617},
            {ReluctanceModels::EFFECTIVE_AREA, 0.017320435},
            {ReluctanceModels::EFFECTIVE_LENGTH, 0.017320435},
            {ReluctanceModels::PARTRIDGE, 0.014844264},
            {ReluctanceModels::STENGLEIN, 0.02528},
            {ReluctanceModels::BALAKRISHNAN, 0.016295213},
            {ReluctanceModels::CLASSIC, 0.02528},
        }));
        INFO("Model: " << magic_enum::enum_name(modelName));
        run_test_energy(modelName, "PQ 40/40", OpenMagneticsTesting::get_ground_gap(0.002), expectedEnergy);
    }

    TEST_CASE("Test_Reluctance_Energy_E_80_38_20_Grinded",
              "[physical-model][reluctance][magnetic-energy]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto [modelName, expectedEnergy] = GENERATE(table<ReluctanceModels, double>({
            {ReluctanceModels::ZHANG, 0.050776555},
            {ReluctanceModels::MUEHLETHALER, 0.04375471},
            {ReluctanceModels::EFFECTIVE_AREA, 0.059827186},
            {ReluctanceModels::EFFECTIVE_LENGTH, 0.059827186},
            {ReluctanceModels::PARTRIDGE, 0.051585101},
            {ReluctanceModels::STENGLEIN, 0.09},
            {ReluctanceModels::BALAKRISHNAN, 0.058340737},
            {ReluctanceModels::CLASSIC, 0.09},
        }));
        INFO("Model: " << magic_enum::enum_name(modelName));
        run_test_energy(modelName, "E 80/38/20", OpenMagneticsTesting::get_ground_gap(0.003), expectedEnergy);
    }

    TEST_CASE("Test_Reluctance_Energy_E_80_38_20_2_Stacks_Grinded",
              "[physical-model][reluctance][magnetic-energy]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto [modelName, expectedEnergy] = GENERATE(table<ReluctanceModels, double>({
            {ReluctanceModels::ZHANG, 0.11118519},
            {ReluctanceModels::MUEHLETHALER, 0.10000631},
            {ReluctanceModels::EFFECTIVE_AREA, 0.12768528},
            {ReluctanceModels::EFFECTIVE_LENGTH, 0.12768528},
            {ReluctanceModels::PARTRIDGE, 0.11476894},
            {ReluctanceModels::STENGLEIN, 0.1629},
            {ReluctanceModels::BALAKRISHNAN, 0.11668147},
            {ReluctanceModels::CLASSIC, 0.1629},
        }));
        INFO("Model: " << magic_enum::enum_name(modelName));
        run_test_energy(modelName, "E 80/38/20", OpenMagneticsTesting::get_ground_gap(0.003), expectedEnergy, 2);
    }

    TEST_CASE("Test_Reluctance_Energy_E_80_38_20_3_Stacks_Grinded",
              "[physical-model][reluctance][magnetic-energy]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto [modelName, expectedEnergy] = GENERATE(table<ReluctanceModels, double>({
            {ReluctanceModels::ZHANG, 0.17222441},
            {ReluctanceModels::MUEHLETHALER, 0.157511},
            {ReluctanceModels::EFFECTIVE_AREA, 0.19591177},
            {ReluctanceModels::EFFECTIVE_LENGTH, 0.19591177},
            {ReluctanceModels::PARTRIDGE, 0.18118045},
            {ReluctanceModels::STENGLEIN, 0.24},
            {ReluctanceModels::BALAKRISHNAN, 0.17502221},
            {ReluctanceModels::CLASSIC, 0.24},
        }));
        INFO("Model: " << magic_enum::enum_name(modelName));
        run_test_energy(modelName, "E 80/38/20", OpenMagneticsTesting::get_ground_gap(0.003), expectedEnergy, 3);
    }

    TEST_CASE("Test_Reluctance_Energy_PQ_40_40_Distributed",
              "[physical-model][reluctance][magnetic-energy]" ALL_RELUCTANCE_MODEL_TAGS "[smoke-test]") {
        auto [modelName, expectedEnergy] = GENERATE(table<ReluctanceModels, double>({
            {ReluctanceModels::ZHANG, 0.045285412},
            {ReluctanceModels::MUEHLETHALER, 0.04483808},
            {ReluctanceModels::EFFECTIVE_AREA, 0.051834892},
            {ReluctanceModels::EFFECTIVE_LENGTH, 0.051834892},
            {ReluctanceModels::PARTRIDGE, 0.046890764},
            {ReluctanceModels::STENGLEIN, 0.031467479},
            {ReluctanceModels::BALAKRISHNAN, 0.050874521},
            {ReluctanceModels::CLASSIC, 0.076},
        }));
        INFO("Model: " << magic_enum::enum_name(modelName));
        run_test_energy(modelName, "PQ 40/40", OpenMagneticsTesting::get_distributed_gap(0.002, 3), expectedEnergy, 1);
    }

    TEST_CASE("Test_Gap_By_Fringing_Factor", "[physical-model][reluctance][gap][smoke-test]") {
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

    TEST_CASE("Test_Reluctance_3C96", "[physical-model][reluctance][smoke-test]") {
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

    TEST_CASE("Test_Web_0", "[physical-model][reluctance][gap][bug][smoke-test]") {
        std::string coreGapData = R"({"area":0.000123,"coordinates":[0,0.0005,0],"distanceClosestNormalSurface":0.014098,"distanceClosestParallelSurface":0.0088,"length":0.000005,"sectionDimensions":[0.0125,0.0125],"shape":"round","type":"subtractive"})";
        std::string modelNameString = "Zhang";
        std::string modelNameStringUpper = modelNameString;
        std::transform(modelNameStringUpper.begin(), modelNameStringUpper.end(), modelNameStringUpper.begin(), ::toupper);
        auto modelName = magic_enum::enum_cast<ReluctanceModels>(modelNameStringUpper);

        auto reluctanceModel = ReluctanceModel::factory(modelName.value());

        CoreGap coreGap(json::parse(coreGapData));

        auto coreGapResult = reluctanceModel->get_gap_reluctance(coreGap);
        // Repro point: the gap reluctance for this web gap must be finite and positive,
        // with a physical fringing factor (>= 1).
        CHECK(std::isfinite(coreGapResult.get_reluctance()));
        CHECK(coreGapResult.get_reluctance() > 0);
        CHECK(coreGapResult.get_fringing_factor() >= 1);
    }

}  // namespace
