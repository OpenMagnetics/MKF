#include "constructive_models/Core.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;


SUITE(CoreProcessedDescription) {
    auto masPath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("MAS/").string();
    auto settings = Settings::GetInstance();

    TEST(E_55_21) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        Core core(coreJson, true);

        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_55_21_N97_additive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000353 * numberStacks, 0.000353 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.124,
                    0.124 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    4.4e-05 * numberStacks, 4.4e-05 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.00035 * numberStacks, 0.00035 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.037, 0.037 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.01015, 0.01015 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0172, 0.0172 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.021 * numberStacks,
                    0.021 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.00935, 0.00935 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.021 * numberStacks,
                    0.021 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.00935, 0.00935 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.021 * numberStacks,
                    0.021 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(E_55_28_21) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_28_21_3C95_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_55_28_21_3C95_additive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000353 * numberStacks, 0.000353 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.124,
                    0.124 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    4.4e-05 * numberStacks, 4.4e-05 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.00035 * numberStacks, 0.00035 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.037, 0.037 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.01015, 0.01015 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0172, 0.0172 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.021 * numberStacks,
                    0.021 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.00935, 0.00935 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.021 * numberStacks,
                    0.021 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.00935, 0.00935 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.021 * numberStacks,
                    0.021 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(E_19_8_5) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_19_8_5_N87_substractive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_19_8_5_N87_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000225 * numberStacks, 0.0000225 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0396,
                    0.0396 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000891 * numberStacks, 0.000000891 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000221 * numberStacks, 0.0000221 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0114, 0.0114 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.00475, 0.00475 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0048, 0.0048 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0048 * numberStacks,
                    0.0048 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.00235, 0.00235 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0048 * numberStacks,
                    0.0048 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.00235, 0.00235 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0048 * numberStacks,
                    0.0048 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(ETD_39_20_13) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_ETD_39_20_13_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ETD 39/20/13";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_ETD_39_20_13_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000125 * numberStacks, 0.000125 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0922,
                    0.0922 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000011500 * numberStacks, 0.000011500 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000123 * numberStacks, 0.000123 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0282, 0.0282 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.00825, 0.00825 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0128, 0.0128 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0128 * numberStacks,
                    0.0128 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0048, 0.0048 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0128 * numberStacks,
                    0.0128 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.0048, 0.0048 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0128 * numberStacks,
                    0.0128 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(ETD_19_14_8) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_ETD_19_14_8_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ETD 19/14/8";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_ETD_19_14_8_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000441 * numberStacks, 0.0000441 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0553,
                    0.0553 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000002440 * numberStacks, 0.000002440 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000395 * numberStacks, 0.0000395 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0184, 0.0184 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0034, 0.0034 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0076, 0.0076 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0076 * numberStacks,
                    0.0076 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.00255, 0.00255 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0076 * numberStacks,
                    0.0076 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.00255, 0.00255 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0076 * numberStacks,
                    0.0076 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(ETD_54_28_19) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_ETD_54_28_19_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ETD 54/28/19";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_ETD_54_28_19_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000280 * numberStacks, 0.000280 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.127,
                    0.127 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000035600 * numberStacks, 0.000035600 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000280 * numberStacks, 0.000280 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0396, 0.0396 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0104, 0.0104 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0193, 0.0193 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0193 * numberStacks,
                    0.0193 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0072, 0.0072 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0193 * numberStacks,
                    0.0193 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.0072, 0.0072 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0193 * numberStacks,
                    0.0193 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(ER_54_18_18) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_ER_54_18_18_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ER 54/18/18";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_ER_54_18_18_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000256 * numberStacks, 0.000256 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.090,
                    0.090 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000023000 * numberStacks, 0.000023000 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000252 * numberStacks, 0.000252 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0216, 0.0216 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.01025, 0.01025 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0183, 0.0183 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0183 * numberStacks,
                    0.0183 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.007, 0.007 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0183 * numberStacks,
                    0.0183 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.007, 0.007 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0183 * numberStacks,
                    0.0183 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(ER_18_3_10) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_ER_18_3_10_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ER 18/3/10";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_ER_18_3_10_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000302 * numberStacks, 0.0000302 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0221,
                    0.0221 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000667 * numberStacks, 0.000000667 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000301 * numberStacks, 0.0000301 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0031, 0.0031 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0047, 0.0047 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0062, 0.0062 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0062 * numberStacks,
                    0.0062 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_minimum_width().value(), 0.0012, 0.0012 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.01 * numberStacks,
                    0.01 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_minimum_width().value(), 0.0012, 0.0012 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.01 * numberStacks,
                    0.01 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(E_102_20_38) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_E_102_20_38_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ELP 102/20/38";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 2;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_102_20_38_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000538 * numberStacks, 0.000538 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.1476,
                    0.1476 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000079410 * numberStacks, 0.000079410 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0005245 * numberStacks, 0.0005245 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0266, 0.0266 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.036, 0.036 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.014, 0.014 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0375 * numberStacks,
                    0.0375 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.008, 0.008 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0375 * numberStacks,
                    0.0375 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.008, 0.008 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0375 * numberStacks,
                    0.0375 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(E_14_3_5_5) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_E_14_3.5_5_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "ELP 14/3.5/5";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 3;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_14_3.5_5_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000143 * numberStacks, 0.0000143 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0207,
                    0.0207 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000296 * numberStacks, 0.000000296 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000139 * numberStacks, 0.0000139 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.004, 0.004 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.004, 0.004 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.003, 0.003 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.005 * numberStacks,
                    0.005 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0015, 0.0015 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.005 * numberStacks,
                    0.005 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.0015, 0.0015 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.005 * numberStacks,
                    0.005 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(EL_25_4_3) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_E_25_4.3_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EL 25/4.3";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_25_4.3_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000856 * numberStacks, 0.0000856 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.030,
                    0.030 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000002570 * numberStacks, 0.000002570 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000083 * numberStacks, 0.000083 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.004, 0.004 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.007255,
                    0.007255 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.00632, 0.00632 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.01454 * numberStacks,
                    0.01454 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.002085, 0.002085 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.002085, 0.002085 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(EL_11_2) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_E_11_2_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EL 11/2";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_E_11_2_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000165 * numberStacks, 0.0000165 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0137,
                    0.0137 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000226 * numberStacks, 0.000000226 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000159 * numberStacks, 0.0000159 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.002, 0.002 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.003195,
                    0.003195 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.00278, 0.00278 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0064 * numberStacks,
                    0.0064 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.000915, 0.000915 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0088 * numberStacks,
                    0.0088 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.000915, 0.000915 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0088 * numberStacks,
                    0.0088 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(EC_70) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EC_70_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EC 70";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EC_70_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000280 * numberStacks, 0.000280 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.144,
                    0.144 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000040420 * numberStacks, 0.000040420 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000211 * numberStacks, 0.000211 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0455, 0.0455 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0141, 0.0141 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0164, 0.0164 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0164 * numberStacks,
                    0.0164 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.013, 0.013 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0164 * numberStacks,
                    0.0164 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.013, 0.013 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0164 * numberStacks,
                    0.0164 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(EFD_10_5_3) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EFD_10_5_3_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EFD 10/5/3";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EFD_10_5_3_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000072 * numberStacks, 0.0000072 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0231,
                    0.0231 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000166 * numberStacks, 0.000000166 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000065 * numberStacks, 0.0000065 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0075, 0.0075 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.00155, 0.00155 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.00455, 0.00455 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.00145 * numberStacks,
                    0.00145 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.001425, 0.001425 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0027 * numberStacks,
                    0.0027 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.001425, 0.001425 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0027 * numberStacks,
                    0.0027 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(EFD_30_15_9) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EFD_30_15_9_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EFD 30/15/9";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EFD_30_15_9_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000069 * numberStacks, 0.000069 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.068,
                    0.068 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000004690 * numberStacks, 0.000004690 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000069 * numberStacks, 0.000069 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0224, 0.0224 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0039, 0.0039 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0146, 0.0146 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0049 * numberStacks,
                    0.0049 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0038, 0.0038 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0091 * numberStacks,
                    0.0091 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.0038, 0.0038 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.0091 * numberStacks,
                    0.0091 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(EQ_30_8_20) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EQ_30_8_20_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EQ 30/8/20";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EQ_30_8_20_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000108 * numberStacks, 0.000108 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.046,
                    0.046 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000004970 * numberStacks, 0.000004970 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000095 * numberStacks, 0.000095 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0106, 0.0106 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0075, 0.0075 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.011, 0.011 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.011 * numberStacks,
                    0.011 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_minimum_width().value(), 0.002, 0.002 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_minimum_width().value(), 0.002, 0.002 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(EPX_10) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EPX_10_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EPX 10";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EPX_10_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000159 * numberStacks, 0.0000159 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0217,
                    0.0217 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000345 * numberStacks, 0.000000345 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000132 * numberStacks, 0.0000132 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0072, 0.0072 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.002825,
                    0.002825 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 3u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.00345, 0.00345 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.00485 * numberStacks,
                    0.00485 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.023 * numberStacks,
                    0.023 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.023 * numberStacks,
                    0.023 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(EPX_7) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EPX_7_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EPX 7";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EPX_7_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000172 * numberStacks, 0.0000172 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0157,
                    0.0157 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000270 * numberStacks, 0.000000270 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000139 * numberStacks, 0.0000139 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0045, 0.0045 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0019, 0.0019 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0034, 0.0034 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0057 * numberStacks,
                    0.0057 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0014 * numberStacks,
                    0.0014 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(EPO_13) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EPO_13_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EPO 13";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EPO_13_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000193 * numberStacks, 0.0000193 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0258,
                    0.0258 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000498 * numberStacks, 0.000000498 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000149 * numberStacks, 0.0000149 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.009, 0.009 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0026, 0.0026 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 3u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0045, 0.0045 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0045 * numberStacks,
                    0.0045 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.026 * numberStacks,
                    0.026 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.026 * numberStacks,
                    0.026 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(LP_42_25_15_8) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_LP_42_25_15_8_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "LP 42/25/15.8";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_LP_42_25_15_8_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000206 * numberStacks, 0.000206 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0901,
                    0.0901 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000018560 * numberStacks, 0.000018560 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000206 * numberStacks, 0.000206 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0228, 0.0228 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0095, 0.0095 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 3u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0162, 0.0162 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0162 * numberStacks,
                    0.0162 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.108 * numberStacks,
                    0.108 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.108 * numberStacks,
                    0.108 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(EP_7) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EP_7_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EP 7";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EP_7_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000103 * numberStacks, 0.0000103 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0157,
                    0.0157 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000162 * numberStacks, 0.000000162 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000085 * numberStacks, 0.0000085 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.005, 0.005 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0019, 0.0019 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0034, 0.0034 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0034 * numberStacks,
                    0.0034 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.048, 0.048 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.00095 * numberStacks,
                    0.00095 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(EP_20) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_EP_20_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "EP 20";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_EP_20_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000078 * numberStacks, 0.000078 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.040,
                    0.040 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000003120 * numberStacks, 0.000003120 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000060 * numberStacks, 0.000060 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.014, 0.014 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.00355, 0.00355 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.009, 0.009 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.009 * numberStacks,
                    0.009 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.129, 0.129 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.00275 * numberStacks,
                    0.00275 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(RM_14) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_RM_14_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "RM 14";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_RM_14_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000200 * numberStacks, 0.000200 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.070,
                    0.070 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000014000 * numberStacks, 0.000014000 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000170 * numberStacks, 0.000170 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0208, 0.0208 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.007, 0.007 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.015, 0.015 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.015 * numberStacks,
                    0.015 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0066 * numberStacks,
                    0.0066 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.0066 * numberStacks,
                    0.0066 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(RM_7LP) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_RM_7LP_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "RM 7LP";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_RM_7LP_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000040 * numberStacks, 0.000040 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0235,
                    0.0235 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000001190 * numberStacks, 0.000001190 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000323 * numberStacks, 0.0000323 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0047, 0.0047 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.00375, 0.00375 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.00725, 0.00725 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.00725 * numberStacks,
                    0.00725 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.002775 * numberStacks,
                    0.002775 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.002775 * numberStacks,
                    0.002775 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(PQ_20_16) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_PQ_20_16_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "PQ 20/16";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_PQ_20_16_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000632 * numberStacks, 0.0000632 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0372,
                    0.0372 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000002360 * numberStacks, 0.000002360 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000544 * numberStacks, 0.0000544 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.01030, 0.01030 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0046, 0.0046 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0088, 0.0088 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0088 * numberStacks,
                    0.0088 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.014 * numberStacks,
                    0.014 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.014 * numberStacks,
                    0.014 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(PQ_107_87) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_PQ_107_87_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "PQ 107/87";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_PQ_107_87_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.001428 * numberStacks, 0.001428 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.204,
                    0.204 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000290600 * numberStacks, 0.000290600 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.001320 * numberStacks, 0.001320 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.056, 0.056 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.02635, 0.02635 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.041, 0.041 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.041 * numberStacks,
                    0.041 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.070 * numberStacks,
                    0.070 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_depth(), 0.070 * numberStacks,
                    0.070 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(PM_114_93) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_PM_114_93_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "PM 114/93";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_PM_114_93_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.001720 * numberStacks, 0.001720 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.200,
                    0.200 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000344000 * numberStacks, 0.000344000 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.001380 * numberStacks, 0.001380 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.063, 0.063 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0225, 0.0225 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.043, 0.043 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.043 * numberStacks,
                    0.043 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.013 * numberStacks,
                    0.013 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.013 * numberStacks,
                    0.013 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(P_150_30) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_P_150_30_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "P 150/30";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_P_150_30_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.003580 * numberStacks, 0.003580 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.160,
                    0.160 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.00056600 * numberStacks, 0.00056600 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.002800 * numberStacks, 0.002800 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.030, 0.030 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0325, 0.0325 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.065, 0.065 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.065 * numberStacks,
                    0.065 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.010 * numberStacks,
                    0.010 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.010 * numberStacks,
                    0.00515 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(P_11_7) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_P_11_7_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "P 11/7";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_P_11_7_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.0000162 * numberStacks, 0.0000162 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.0155,
                    0.0155 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000251 * numberStacks, 0.000000251 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.0000132 * numberStacks, 0.0000132 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.00440, 0.00440 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.00215, 0.00215 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0047, 0.0047 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0047 * numberStacks,
                    0.0047 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.00105 * numberStacks,
                    0.00105 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.00105 * numberStacks,
                    0.00105 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(P_7_4) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_P_7_4_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "P 7/4";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_P_7_4_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000007 * numberStacks, 0.000007 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.010,
                    0.010 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000070 * numberStacks, 0.000000070 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000006 * numberStacks, 0.000006 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0028, 0.0028 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0014, 0.0014 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.003, 0.003 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.003 * numberStacks,
                    0.003 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.00075 * numberStacks,
                    0.00075 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[2].get_width(), 0.00075 * numberStacks,
                    0.00075 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
        CHECK(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
    }

    TEST(U_79_129_31) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_U_79_129_31_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "U 79/129/31";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_U_79_129_31_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000693 * numberStacks, 0.000693 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.309,
                    0.309 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000214220 * numberStacks, 0.000214220 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000693 * numberStacks, 0.000693 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.085, 0.085 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.034, 0.034 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.022, 0.022 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0315 * numberStacks,
                    0.0315 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.022 * numberStacks,
                    0.022 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0315 * numberStacks,
                    0.0315 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(U_26_22_16) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_U_26_22_16_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "U 26/22/16";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_U_26_22_16_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000131 * numberStacks, 0.000131 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.098,
                    0.098 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000012800 * numberStacks, 0.000012800 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000129 * numberStacks, 0.000129 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.026, 0.026 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.009, 0.009 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0084, 0.0084 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.016 * numberStacks,
                    0.016 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0084 * numberStacks,
                    0.0084 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.016 * numberStacks,
                    0.016 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(UR_48_39_17) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_UR_48_39_17_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "UR 48/39/17";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_UR_48_39_17_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000215 * numberStacks, 0.000215 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.186,
                    0.186 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000039990 * numberStacks, 0.000039990 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000215 * numberStacks, 0.000215 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0538, 0.0538 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0174, 0.0174 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.017, 0.017 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.017 * numberStacks,
                    0.017 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.013 * numberStacks,
                    0.013 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.017 * numberStacks,
                    0.017 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(UR_70_33_17) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_UR_70_33_17_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "UR 70/33/17";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_UR_70_33_17_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000214 * numberStacks, 0.000214 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.197,
                    0.197 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000043800 * numberStacks, 0.000043800 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000214 * numberStacks, 0.000214 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.0381, 0.0381 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.035, 0.035 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.01725, 0.01725 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.01725 * numberStacks,
                    0.01725 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.01725 * numberStacks,
                    0.01725 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.01725 * numberStacks,
                    0.01725 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::ROUND);
    }

    TEST(UR_55_39_36) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_UR_55_39_36_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "UR 55/38/36";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_UR_55_39_36_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000418 * numberStacks, 0.000418 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.188,
                    0.188 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000078570 * numberStacks, 0.000078570 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000418 * numberStacks, 0.000418 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.051, 0.051 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0196, 0.0196 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0235, 0.0235 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0235 * numberStacks,
                    0.0235 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.012 * numberStacks,
                    0.012 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.036 * numberStacks,
                    0.036 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(UR_64_40_20) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_UR_64_40_20_N97_substractive";
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "UR 64/40/20";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_UR_64_40_20_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000290 * numberStacks, 0.000290 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.210,
                    0.210 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000061000 * numberStacks, 0.000061000 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000290 * numberStacks, 0.000290 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.053, 0.053 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0232, 0.0232 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.020, 0.020 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.020 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::ROUND);
    }

    TEST(UT_20) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_UT_20_N97";
        coreJson["functionalDescription"]["type"] = "closed shape";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "UT 20";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        auto geometrical_description = *(core.get_geometrical_description());

        CHECK_EQUAL(*(core.get_name()), "core_UT_20_N97");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::CLOSED);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000013 * numberStacks, 0.000013 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.053,
                    0.053 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000000688 * numberStacks, 0.000000688 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000013 * numberStacks, 0.000013 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_height()), 0.016, 0.016 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_width()), 0.0075, 0.0075 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.0041, 0.0041 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.0046 * numberStacks,
                    0.0046 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_width(), 0.0033 * numberStacks,
                    0.0033 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[1].get_depth(), 0.0046 * numberStacks,
                    0.0046 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() ==
              ColumnShape::RECTANGULAR);
        CHECK(core.get_processed_description()->get_columns()[1].get_shape() ==
              ColumnShape::RECTANGULAR);
    }

    TEST(T_40_24_16) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_T_40_24_16_N97_substractive";
        coreJson["functionalDescription"]["type"] = "toroidal";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "T 40/24/16";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        Core core(coreJson, true);
        double numberStacks = coreJson["functionalDescription"]["numberStacks"];

        CHECK_EQUAL(*(core.get_name()), "core_T_40_24_16_N97_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(),
                    0.000125 * numberStacks, 0.000125 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_length(), 0.09629,
                    0.09629 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_volume(),
                    0.000012070 * numberStacks, 0.000012070 * numberStacks * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_minimum_area(),
                    0.000125 * numberStacks, 0.000125 * numberStacks * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_radial_height()), 0.012, 0.012 * 0.2);
        CHECK_CLOSE(*(core.get_processed_description()->get_winding_windows()[0].get_angle()), 360, 360 * 0.2);
        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 1u);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_width(), 0.008, 0.008 * 0.2);
        CHECK_CLOSE(core.get_processed_description()->get_columns()[0].get_depth(), 0.016 * numberStacks,
                    0.020 * numberStacks * 0.2);
        CHECK(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::RECTANGULAR);
    }

    TEST(Core_Procecessd_Description_Web_0) {
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 0.000123, \"coordinates\": "
            "[0.0, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.01455, \"length\": 0.0001, \"sectionDimensions\": "
            "[0.0125, 0.0125], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 6.2e-05, \"coordinates\": "
            "[0.017301, 0.0005, 0.0], \"distanceClosestNormalSurface\": 0.014598, \"length\": 5e-06, "
            "\"sectionDimensions\": [0.004501, 0.0125], \"shape\": \"irregular\", \"type\": \"residual\"}, {\"area\": "
            "6.2e-05, \"coordinates\": [-0.017301, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.014598, \"length\": "
            "5e-06, \"sectionDimensions\": [0.004501, 0.0125], \"shape\": \"irregular\", \"type\": \"residual\"}], "
            "\"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.0125, "
            "\"B\": 0.0064, \"C\": 0.0088, \"D\": 0.0046, \"E\": 0.01, \"F\": 0.0043, \"G\": 0.0, \"H\": 0.0, \"K\": "
            "0.0023}, \"family\": \"ep\", \"familySubtype\": \"1\", \"name\": \"Custom\", \"type\": \"custom\"}, "
            "\"type\": \"two-piece set\"}, \"geometricalDescription\": null, \"processedDescription\": null}");

        Core core(coreJson, true);

        CHECK_EQUAL(core.get_processed_description()->get_columns().size(), 2u);
    }

    TEST(Core_Procecessd_Description_Web_1) {
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [], \"material\": \"3C97\", "
            "\"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.0308, \"B\": 0.0264, \"C\": "
            "0.0265, \"D\": 0.016, \"E\": 0.01, \"G\": 0.0, \"H\": 0.0}, \"family\": \"u\", \"familySubtype\": \"1\", "
            "\"name\": \"Custom\", \"type\": \"custom\"}, \"type\": \"two-piece set\"}}");

        Core core(coreJson, true);

        CHECK_CLOSE(core.get_processed_description()->get_effective_parameters().get_effective_area(), 0.0002756,
                    0.0002756 * 0.2);
        auto functionalDescription = core.get_functional_description();
    }


    TEST(Test_Core_All_Shapes) {
        settings->set_use_toroidal_cores(true);
        auto shapeNames = get_shape_names();
        for (auto shapeName : shapeNames) {
            if (shapeName.contains("PQI") || shapeName.contains("UI ")) {
                continue;
            }
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
            if ((core.get_processed_description()->get_effective_parameters().get_effective_area() <= 0) ||
                (core.get_processed_description()->get_effective_parameters().get_effective_length() <= 0) ||
                (core.get_processed_description()->get_effective_parameters().get_effective_volume() <= 0) ||
                (core.get_processed_description()->get_effective_parameters().get_minimum_area() <= 0)) {

                std::cout << "shapeName: " << shapeName << std::endl;
            }

            CHECK(core.get_processed_description()->get_effective_parameters().get_effective_area() > 0);
            CHECK(core.get_processed_description()->get_effective_parameters().get_effective_length() > 0);
            CHECK(core.get_processed_description()->get_effective_parameters().get_effective_volume() > 0);
            CHECK(core.get_processed_description()->get_effective_parameters().get_minimum_area() > 0);
        }
    }
}

SUITE(CoreGeometricalDescription) {
    auto settings = Settings::GetInstance();
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");

    TEST(E_19_8_5_Geometrical_Description) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_19_8_5_N87_substractive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        Core core(coreJson, true);

        auto geometrical_description = *(core.get_geometrical_description());

        CHECK_EQUAL(*(core.get_name()), "core_E_19_8_5_N87_substractive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_EQUAL(geometrical_description.size(), 2u);
        CHECK(geometrical_description[0].get_machining());
        CHECK(!geometrical_description[1].get_machining());
        CHECK(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
        CHECK(geometrical_description[1].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
    }

    TEST(E_55_21) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        Core core(coreJson, true);

        auto geometrical_description = *(core.get_geometrical_description());

        CHECK_EQUAL(*(core.get_name()), "core_E_55_21_N97_additive");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_EQUAL(geometrical_description.size(), 6u);
        CHECK(!geometrical_description[0].get_machining());
        CHECK(!geometrical_description[1].get_machining());
        CHECK(!geometrical_description[2].get_machining());
        CHECK(!geometrical_description[3].get_machining());
        CHECK(!geometrical_description[4].get_machining());
        CHECK(!geometrical_description[5].get_machining());
        CHECK(geometrical_description[4].get_type() == CoreGeometricalDescriptionElementType::SPACER);
        CHECK(geometrical_description[5].get_type() == CoreGeometricalDescriptionElementType::SPACER);
    }

    TEST(T_40_24_16) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_T_40_24_16_N97.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        Core core(coreJson, true);

        auto geometrical_description = *(core.get_geometrical_description());

        CHECK_EQUAL(*(core.get_name()), "core_T_40_24_16_N97");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_EQUAL(geometrical_description.size(), 1u);
    }

    TEST(Core_Web_0) {
        auto coreJson = json::parse(R"(
            {"name": "Custom_0", "functionalDescription": {"gapping": [{"area": 0.000114, "coordinates":
            [0.0, -0.00425, 0.0], "distanceClosestNormalSurface": 0.004201, "length": 0.0001,
            "sectionDimensions": [0.012, 0.012], "shape": "round", "type": "subtractive"}, {"area":
            0.000114, "coordinates": [0.0, 0.0, 0.0], "distanceClosestNormalSurface": 0.008451, "length":
            0.0001, "sectionDimensions": [0.012, 0.012], "shape": "round", "type": "subtractive"},
            {"area": 0.000114, "coordinates": [0.0, 0.00425, 0.0], "distanceClosestNormalSurface": 0.004201,
            "length": 0.0001, "sectionDimensions": [0.012, 0.012], "shape": "round", "type":
            "subtractive"}, {"area": 0.000205, "coordinates": [0.017925, 0.0, 0.0],
            "distanceClosestNormalSurface": 0.0085, "length": 5e-06, "sectionDimensions": [0.01025, 0.02],
            "shape": "irregular", "type": "residual"}, {"area": 0.000205, "coordinates": [-0.017925, 0.0,
            0.0], "distanceClosestNormalSurface": 0.0085, "length": 5e-06, "sectionDimensions": [0.01025,
            0.02], "shape": "irregular", "type": "residual"}], "material": "3C97", "numberStacks": 1,
            "shape": {"aliases": [], "dimensions": {"A": 0.03, "B": 0.011800000000000001, "C": 0.02,
            "D": 0.0085, "E": 0.0256, "F": 0.012, "G": 0.017, "H": 0.0}, "family": "lp",
            "familySubtype": "1", "name": "Custom", "type": "custom"}, "type": "two-piece set"}})");

        Core core(coreJson, true);

        auto geometrical_description = *(core.get_geometrical_description());

        CHECK_EQUAL(*(core.get_name()), "Custom_0");
        CHECK(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
                  .get_mutable_volumetric_losses()["default"]
                  .size() > 0);
        CHECK_EQUAL(geometrical_description.size(), 2u);
        CHECK(geometrical_description[0].get_machining());
        CHECK(geometrical_description[1].get_machining());
        CHECK(geometrical_description[0].get_machining()->size() == 2);
        CHECK(geometrical_description[1].get_machining()->size() == 2);
        CHECK(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
        CHECK(geometrical_description[1].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
    }

    TEST(Test_Core_Geometrical_Description_Web_1) {
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 0.000135, \"coordinates\": "
            "[0.0, 0.0078, 0.0], \"distanceClosestNormalSurface\": 0.00515, \"length\": 0.0001, \"sectionDimensions\": "
            "[0.008401, 0.016], \"shape\": \"rectangular\", \"type\": \"subtractive\"}, {\"area\": 0.000135, "
            "\"coordinates\": [0.0, 0.0026, 0.0], \"distanceClosestNormalSurface\": 0.0047, \"length\": 0.001, "
            "\"sectionDimensions\": [0.008401, 0.016], \"shape\": \"rectangular\", \"type\": \"subtractive\"}, "
            "{\"area\": 0.000135, \"coordinates\": [0.0, -0.0020299999999999997, 0.0], "
            "\"distanceClosestNormalSurface\": 0.00512, \"length\": 0.00016, \"sectionDimensions\": [0.008401, 0.016], "
            "\"shape\": \"rectangular\", \"type\": \"subtractive\"}, {\"area\": 0.000135, \"coordinates\": [0.0, "
            "-0.007549999999999999, 0.0], \"distanceClosestNormalSurface\": 0.0027, \"length\": 0.005, "
            "\"sectionDimensions\": [0.008401, 0.016], \"shape\": \"rectangular\", \"type\": \"subtractive\"}, "
            "{\"area\": 0.000135, \"coordinates\": [0.0174, 0.005, 0.0], \"distanceClosestNormalSurface\": 0.008, "
            "\"length\": 0.003, \"sectionDimensions\": [0.008401, 0.016], \"shape\": \"rectangular\", \"type\": "
            "\"subtractive\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], "
            "\"dimensions\": {\"A\": 0.0258, \"B\": 0.0222, \"C\": 0.016, \"D\": 0.013, \"E\": 0.009, \"F\": 0.0125, "
            "\"G\": 0.0, \"H\": 0.0}, \"family\": \"u\", \"familySubtype\": \"1\", \"name\": \"Custom\", \"type\": "
            "\"custom\"}, \"type\": \"two-piece set\"}, \"geometricalDescription\": null, \"processedDescription\": "
            "null}");

        Core core(coreJson, true);

        auto geometrical_description = *(core.get_geometrical_description());

        CHECK_EQUAL(geometrical_description.size(), 2u);
        CHECK(geometrical_description[0].get_machining());
        CHECK(geometrical_description[1].get_machining());
        CHECK(geometrical_description[0].get_machining()->size() == 3);
        CHECK(geometrical_description[1].get_machining()->size() == 2);
        CHECK(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
        CHECK(geometrical_description[1].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
    }
}

SUITE(CoreFunctionalDescription) {
    auto settings = Settings::GetInstance();
    std::string filePath = __FILE__;
    double maximumError = 0.05;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");

    TEST(E_55_21_all_gaps_residual) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        coreJson["functionalDescription"]["gapping"][0] = coreJson["functionalDescription"]["gapping"][1];

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 3u);
        CHECK(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
        CHECK(*(functionalDescription.get_gapping()[0].get_shape()) ==
              *(functionalDescription.get_gapping()[1].get_shape()));
        CHECK(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) ==
              *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
        CHECK(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[1].get_area()) * 2, 0.2);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[0] ==
              -(*functionalDescription.get_gapping()[2].get_coordinates())[0]);
    }

    TEST(E_55_21_central_gap) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 3u);
        CHECK(functionalDescription.get_gapping()[0].get_type() != functionalDescription.get_gapping()[1].get_type());
        CHECK(*(functionalDescription.get_gapping()[0].get_shape()) ==
              *(functionalDescription.get_gapping()[1].get_shape()));
        CHECK(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) !=
              *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
        CHECK(functionalDescription.get_gapping()[0].get_length() != functionalDescription.get_gapping()[1].get_length());
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[1].get_area()) * 2, 0.2);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] != 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[0] ==
              -(*functionalDescription.get_gapping()[2].get_coordinates())[0]);
    }

    TEST(E_55_21_gap_all_columns) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        coreJson["functionalDescription"]["gapping"][1] = coreJson["functionalDescription"]["gapping"][0];
        coreJson["functionalDescription"]["gapping"][2] = coreJson["functionalDescription"]["gapping"][0];

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 3u);
        CHECK(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
        CHECK(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[2].get_type());
        CHECK(*(functionalDescription.get_gapping()[0].get_shape()) ==
              *(functionalDescription.get_gapping()[1].get_shape()));
        CHECK(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) ==
              *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
        CHECK(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[1].get_area()) * 2, 0.2);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[0] ==
              -(*functionalDescription.get_gapping()[2].get_coordinates())[0]);
    }

    TEST(E_55_21_central_distributed_gap_even) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        coreJson["functionalDescription"]["gapping"].push_back(coreJson["functionalDescription"]["gapping"][0]);

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 4u);
        CHECK(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
        CHECK(functionalDescription.get_gapping()[0].get_type() != functionalDescription.get_gapping()[2].get_type());
        CHECK(*(functionalDescription.get_gapping()[0].get_shape()) ==
              *(functionalDescription.get_gapping()[1].get_shape()));
        CHECK(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) ==
              *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
        CHECK(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) !=
              *(functionalDescription.get_gapping()[2].get_distance_closest_normal_surface()));
        CHECK(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
        CHECK(functionalDescription.get_gapping()[0].get_length() != functionalDescription.get_gapping()[2].get_length());
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[1].get_area()), 0.2);
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[2].get_area()) * 2, 0.2);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] ==
              (*functionalDescription.get_gapping()[1].get_coordinates())[0]);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] ==
              -(*functionalDescription.get_gapping()[1].get_coordinates())[1]);
        CHECK((*functionalDescription.get_gapping()[2].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[2].get_coordinates())[0] ==
              -(*functionalDescription.get_gapping()[3].get_coordinates())[0]);
    }

    TEST(E_55_21_central_distributed_gap_odd) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);
        coreJson["functionalDescription"]["gapping"].push_back(coreJson["functionalDescription"]["gapping"][0]);
        coreJson["functionalDescription"]["gapping"].push_back(coreJson["functionalDescription"]["gapping"][0]);

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 5u);
        CHECK(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
        CHECK(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[2].get_type());
        CHECK(functionalDescription.get_gapping()[0].get_type() != functionalDescription.get_gapping()[3].get_type());
        CHECK(*(functionalDescription.get_gapping()[0].get_shape()) ==
              *(functionalDescription.get_gapping()[1].get_shape()));
        CHECK(*(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()) >
              *(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()));
        CHECK(*(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()) >
              *(functionalDescription.get_gapping()[2].get_distance_closest_normal_surface()));
        CHECK(*(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()) <
              *(functionalDescription.get_gapping()[3].get_distance_closest_normal_surface()));
        CHECK(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
        CHECK(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[2].get_length());
        CHECK(functionalDescription.get_gapping()[0].get_length() != functionalDescription.get_gapping()[3].get_length());
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[1].get_area()), 0.2);
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[2].get_area()), 0.2);
        CHECK_CLOSE(*(functionalDescription.get_gapping()[0].get_area()),
                    *(functionalDescription.get_gapping()[3].get_area()) * 2, 0.2);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] ==
              (*functionalDescription.get_gapping()[1].get_coordinates())[0]);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] ==
              (*functionalDescription.get_gapping()[2].get_coordinates())[0]);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] ==
              -(*functionalDescription.get_gapping()[2].get_coordinates())[1]);
        CHECK((*functionalDescription.get_gapping()[3].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[3].get_coordinates())[0] ==
              -(*functionalDescription.get_gapping()[4].get_coordinates())[0]);
    }

    TEST(Test_Core_Functional_Description_Web_0) {
        auto coreJson = json::parse(
            ""
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [], \"material\": \"3C97\", \"shape\": "
            "{\"family\": \"pm\", \"type\": \"custom\", \"aliases\": [], \"dimensions\": {\"A\": 0.1118, \"B\": "
            "0.046299999999999994, \"C\": 0.045, \"D\": 0.0319, \"E\": 0.08979999999999999, \"F\": 0.0286, \"G\": "
            "0.052, \"H\": 0.0056, \"b\": 0.0058, \"t\": 0.004200000000000001}, \"familySubtype\": \"2\", \"name\": "
            "\"Custom\"}, \"type\": \"two-piece set\", \"numberStacks\": 1}}"
            "");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
    }

    TEST(Test_Core_Functional_Description_Web_1) {
        // Tests that a missaligned gapping get recalculated
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 0.000123, \"coordinates\": "
            "[0.0, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.01455, \"length\": 0.0001, \"sectionDimensions\": "
            "[0.0125, 0.0125], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 6.2e-05, \"coordinates\": "
            "[0.017301, 0.0005, 0.0], \"distanceClosestNormalSurface\": 0.014598, \"length\": 5e-06, "
            "\"sectionDimensions\": [0.004501, 0.0125], \"shape\": \"irregular\", \"type\": \"residual\"}, {\"area\": "
            "6.2e-05, \"coordinates\": [-0.017301, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.014598, \"length\": "
            "5e-06, \"sectionDimensions\": [0.004501, 0.0125], \"shape\": \"irregular\", \"type\": \"residual\"}], "
            "\"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.0125, "
            "\"B\": 0.0064, \"C\": 0.0088, \"D\": 0.0046, \"E\": 0.01, \"F\": 0.0043, \"G\": 0.0, \"H\": 0.0, \"K\": "
            "0.0023}, \"family\": \"ep\", \"familySubtype\": \"1\", \"name\": \"Custom\", \"type\": \"custom\"}, "
            "\"type\": \"two-piece set\"}, \"geometricalDescription\": null, \"processedDescription\": null}");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 2u);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[2] == 0);

        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[2] != 0);
    }

    TEST(Test_Core_Functional_Description_Web_2) {
        // Tests that a distributed but aligned gapping does not get recalculated
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 1.5e-05, \"coordinates\": "
            "[0.0, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.0041, \"length\": 0.001, \"sectionDimensions\": "
            "[0.0043, 0.0043], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 1.5e-05, \"coordinates\": "
            "[0.0, 0.001, 0.0], \"distanceClosestNormalSurface\": 0.0041, \"length\": 0.001, \"sectionDimensions\": "
            "[0.0043, 0.0043], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 8.8e-05, \"coordinates\": "
            "[0.0, 0.0, -0.005751], \"distanceClosestNormalSurface\": 0.004598, \"length\": 5e-06, "
            "\"sectionDimensions\": [0.058628, 0.001501], \"shape\": \"irregular\", \"type\": \"residual\"}, "
            "{\"area\": 8.8e-05, \"coordinates\": [0.0, -0.001, -0.005751], \"distanceClosestNormalSurface\": 0.004598, "
            "\"length\": 5e-06, \"sectionDimensions\": [0.058628, 0.001501], \"shape\": \"irregular\", \"type\": "
            "\"residual\"}], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.0125, \"B\": 0.0064, \"C\": 0.0088, \"D\": 0.0046, \"E\": 0.01, \"F\": 0.0043, \"G\": 0.000, "
            "\"H\": 0.0, \"K\": 0.0023}, \"family\": \"ep\", \"familySubtype\": \"1\", \"name\": \"Custom\", \"type\": "
            "\"custom\"}, \"type\": \"two-piece set\"}, \"geometricalDescription\": null, \"processedDescription\": "
            "null}");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();

        CHECK_EQUAL(functionalDescription.get_gapping().size(), 4u);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[0].get_coordinates())[2] == 0);

        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[1] != 0);
        CHECK((*functionalDescription.get_gapping()[1].get_coordinates())[2] == 0);

        CHECK((*functionalDescription.get_gapping()[2].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[2].get_coordinates())[1] == 0);
        CHECK((*functionalDescription.get_gapping()[2].get_coordinates())[2] != 0);

        CHECK((*functionalDescription.get_gapping()[3].get_coordinates())[0] == 0);
        CHECK((*functionalDescription.get_gapping()[3].get_coordinates())[1] != 0);
        CHECK((*functionalDescription.get_gapping()[3].get_coordinates())[2] != 0);
    }

    TEST(Test_Core_Functional_Description_Web_3) {
        // Check for segmentation fault
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 1.5e-05, \"coordinates\": "
            "[0.0, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.00455, \"length\": 0.0001, \"sectionDimensions\": "
            "[0.0043, 0.0043], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 8.8e-05, \"coordinates\": "
            "[0.0, 0.0, -0.005751], \"distanceClosestNormalSurface\": 0.004598, \"length\": 5e-06, "
            "\"sectionDimensions\": [0.058628, 0.001501], \"shape\": \"irregular\", \"type\": \"residual\"}], "
            "\"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.101, "
            "\"B\": 0.076, \"C\": 0.03, \"D\": 0.048, \"E\": 0.044, \"G\": 0.0, \"H\": 0.0}, \"family\": \"u\", "
            "\"familySubtype\": \"1\", \"name\": \"Custom\", \"type\": \"custom\"}, \"type\": \"two-piece set\"}, "
            "\"geometricalDescription\": null, \"processedDescription\": null}");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
    }

    TEST(Test_Core_Functional_Description_Web_4) {
        // Check for segmentation fault
        auto coreJson = json::parse(
            "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 0.000175, \"coordinates\": "
            "[0.0, -0.0124, 0.0], \"distanceClosestNormalSurface\": 0.0119, \"length\": 0.001, \"sectionDimensions\": "
            "[0.0149, 0.0149], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 0.000175, \"coordinates\": "
            "[0.0, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.024301, \"length\": 0.002, \"sectionDimensions\": "
            "[0.0149, 0.0149], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 0.000175, \"coordinates\": "
            "[0.0, 0.0124, 0.0], \"distanceClosestNormalSurface\": 0.011901, \"length\": 0.002, \"sectionDimensions\": "
            "[0.0149, 0.0149], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 0.000136, \"coordinates\": "
            "[0.0344, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.0248, \"length\": 5e-06, \"sectionDimensions\": "
            "[0.0091, 0.0149], \"shape\": \"rectangular\", \"type\": \"residual\"}], \"material\": \"3C97\", "
            "\"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.038700000000000005, \"B\": "
            "0.0352, \"C\": 0.0149, \"D\": 0.0248, \"G\": 0.0, \"H\": 0.0091}, \"family\": \"ur\", \"familySubtype\": "
            "\"1\", \"name\": \"Custom\", \"type\": \"custom\"}, \"type\": \"two-piece set\"}, "
            "\"geometricalDescription\": null, \"processedDescription\": null}");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
    }

    TEST(Test_Core_Functional_Description_Web_5) {
        // Check for segmentation fault
        auto coreJson = json::parse(
            "{\"name\": \"dummy\", \"functionalDescription\": {\"gapping\": [{\"length\": 0.001, \"type\": "
            "\"subtractive\"}, {\"length\": 0.002, \"type\": \"subtractive\"}, {\"length\": 0.002, \"type\": "
            "\"subtractive\"}, {\"length\": 0.00005, \"type\": \"residual\"}, {\"length\": 0.00005, \"type\": "
            "\"residual\"}], \"material\": \"N97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
            "{\"A\": 0.0112, \"B\": 0.0052, \"C\": 0.0045000000000000005, \"D\": 0.0036, \"E\": 0.008150000000000001, "
            "\"F\": 0.0038, \"G\": 0.0058, \"H\": 0.0020499999999999997, \"J\": 0.009600000000000001, \"R\": 0.0003}, "
            "\"family\": \"rm\", \"familySubtype\": \"3\", \"magneticCircuit\": \"open\", \"name\": \"RM 4\", "
            "\"type\": \"standard\"}, \"type\": \"two-piece set\"}}");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
        CHECK_EQUAL(functionalDescription.get_gapping().size(), 5u);
    }

    TEST(Test_Core_Functional_Description_Web_6) {
        // Check for segmentation fault
        auto coreJson = json::parse(
            "{\"name\": \"My Core test 2\", \"functionalDescription\": {\"gapping\": [{\"area\": 0.000199, "
            "\"coordinates\": [0.0, 0.0005, 0.0], \"distanceClosestNormalSurface\": 0.0064, "
            "\"distanceClosestParallelSurface\": 0.0072499999999999995, \"length\": 0.001, \"sectionDimensions\": "
            "[0.015901, 0.015901], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 0.000123, "
            "\"coordinates\": [0.0165, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.007396, "
            "\"distanceClosestParallelSurface\": 0.0072499999999999995, \"length\": 1e-05, \"sectionDimensions\": "
            "[0.0026, 0.047308], \"shape\": \"irregular\", \"type\": \"residual\"}, {\"area\": 0.000123, "
            "\"coordinates\": [-0.0165, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.007396, "
            "\"distanceClosestParallelSurface\": 0.0072499999999999995, \"length\": 1e-05, \"sectionDimensions\": "
            "[0.0026, 0.047308], \"shape\": \"irregular\", \"type\": \"residual\"}], \"material\": \"N92\", "
            "\"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.0577, \"B\": "
            "0.028399999999999998, \"C\": 0.0155, \"D\": 0.016, \"H\": 0.01590}, \"family\": \"ur\", "
            "\"familySubtype\": \"2\", \"magneticCircuit\": null, \"name\": \"UR 57/28/16\", \"type\": \"standard\"}, "
            "\"type\": \"two-piece set\"}, \"geometricalDescription\": null, \"processedDescription\": null}");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
        CHECK_EQUAL(functionalDescription.get_gapping().size(), 2u);
    }

    TEST(Test_Core_Functional_Description_Web_7) {
        // Check for segmentation fault
        auto coreJson = json::parse(R"({"name":"My Core","functionalDescription":{"coating":null,"gapping":[{"area":0.000057,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.010097499999999999,"distanceClosestParallelSurface":0.005050000000000001,"length":0.000005,"sectionDimensions":[0.0085,0.0085],"shape":"round","type":"residual"},{"area":0.000028,"coordinates":[0.01075,0,0],"distanceClosestNormalSurface":0.010097499999999999,"distanceClosestParallelSurface":0.005050000000000001,"length":0.000005,"sectionDimensions":[0.0029,0.0085],"shape":"irregular","type":"residual"},{"area":0.000028,"coordinates":[-0.01075,0,0],"distanceClosestNormalSurface":0.010097499999999999,"distanceClosestParallelSurface":0.005050000000000001,"length":0.000005,"sectionDimensions":[0.0029,0.0085],"shape":"irregular","type":"residual"}],"material":"3C97","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":0.0576,"B":0.028399999999999998,"C":0.0155,"D":0.016,"H":0.0159,"G":0},"family":"ur","familySubtype":"2","name":"UR 57/28/16","type":"standard"},"type":"two-piece set"}})");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
        CHECK_EQUAL(functionalDescription.get_gapping().size(), 2u);
    }

    TEST(Test_Core_Functional_Description_Web_8) {
        // Check for segmentation fault
        auto coreJson = json::parse(R"({"functionalDescription": {"type": "two-piece set", "material": "3C97", "shape": "U 80/150/30", "gapping": [{"length": 0.003, "type": "additive", "coordinates": [0, 0, 0 ] }, {"length": 0.003, "type": "additive", "coordinates": [0.0595, 0, 0 ] } ], "numberStacks": 1 }, "name": "My Core", "geometricalDescription": null, "processedDescription": null })");

        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
        CHECK_EQUAL(functionalDescription.get_gapping().size(), 2u);
    }

    TEST(Test_Core_Functional_Description_Web_9) {
        Core core = json::parse(R"({"distributorsInfo":[{"cost":1.17,"country":"USA","distributedArea":"International","email":null,"link":"https://www.digikey.com/en/products/detail/ferroxcube/E18-4-10-R-3F36/7041469","name":"Digi-Key","phone":null,"quantity":7063,"reference":"1779-1009-ND","updatedAt":"05/10/2023"},{"cost":null,"country":"UK","distributedArea":"International","email":null,"link":"https://www.shop.gatewaycando.com/magnetics/cores","name":"Gateway","phone":null,"quantity":5,"reference":"E18/4/10/R-3F36","updatedAt":"05/10/2023"}],"functionalDescription":{"coating":null,"gapping":[{"area":0.00004,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.0019975,"distanceClosestParallelSurface":0.005,"length":0.000005,"sectionDimensions":[0.004,0.01],"shape":"rectangular","type":"residual"},{"area":0.000021,"coordinates":[0.008,0,0],"distanceClosestNormalSurface":0.0019975,"distanceClosestParallelSurface":0.005,"length":0.000005,"sectionDimensions":[0.002001,0.01],"shape":"rectangular","type":"residual"},{"area":0.000021,"coordinates":[-0.008,0,0],"distanceClosestNormalSurface":0.0019975,"distanceClosestParallelSurface":0.005,"length":0.000005,"sectionDimensions":[0.002001,0.01],"shape":"rectangular","type":"residual"}],"material":{"bhCycle":null,"coerciveForce":[{"magneticField":32,"magneticFluxDensity":0,"temperature":100},{"magneticField":37,"magneticFluxDensity":0,"temperature":25}],"curieTemperature":230,"density":4750,"family":"3F","heatCapacity":{"excludeMaximum":null,"excludeMinimum":null,"maximum":800,"minimum":700,"nominal":null},"heatConductivity":{"excludeMaximum":null,"excludeMinimum":null,"maximum":5,"minimum":3.5,"nominal":null},"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Ferroxcube","reference":null,"status":null},"material":"ferrite","name":"3F36","permeability":{"amplitude":null,"initial":[{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-40,"tolerance":null,"value":1577},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-30,"tolerance":null,"value":1590},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-20,"tolerance":null,"value":1611},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-10,"tolerance":null,"value":1633},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":0,"tolerance":null,"value":1657},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":10,"tolerance":null,"value":1683},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":20,"tolerance":null,"value":1710},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":30,"tolerance":null,"value":1767},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":40,"tolerance":null,"value":1792},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":50,"tolerance":null,"value":1818},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":60,"tolerance":null,"value":1842},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":70,"tolerance":null,"value":1780},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":80,"tolerance":null,"value":1794},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":90,"tolerance":null,"value":1807},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":100,"tolerance":null,"value":1818},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":110,"tolerance":null,"value":1893},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":120,"tolerance":null,"value":1897},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":130,"tolerance":null,"value":1901},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":140,"tolerance":null,"value":1907},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":150,"tolerance":null,"value":1848},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":160,"tolerance":null,"value":1855},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":170,"tolerance":null,"value":1865},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":180,"tolerance":null,"value":1878}]},"remanence":[{"magneticField":0,"magneticFluxDensity":0.105,"temperature":100},{"magneticField":0,"magneticFluxDensity":0.125,"temperature":25}],"resistivity":[{"temperature":25,"value":12}],"saturation":[{"magneticField":1200,"magneticFluxDensity":0.42,"temperature":100},{"magneticField":1200,"magneticFluxDensity":0.52,"temperature":25}],"type":"commercial","volumetricLosses":{"default":[{"a":null,"b":null,"c":null,"coefficients":{"excessLossesCoefficient":1.82280296e-20,"resistivityFrequencyCoefficient":4.25672064e-28,"resistivityMagneticFluxDensityCoefficient":14.7587264,"resistivityOffset":6.465016450000001e-17,"resistivityTemperatureCoefficient":5.84938089e-16},"d":null,"method":"roshen","ranges":null,"referenceVolumetricLosses":null},{"a":null,"b":null,"c":null,"coefficients":null,"d":null,"method":"steinmetz","ranges":[{"alpha":1.43902,"beta":3.26718,"ct0":1.232717265,"ct1":0.010783518,"ct2":0.00008394600000000001,"k":6.83,"maximumFrequency":499999,"minimumFrequency":100000},{"alpha":2.19515,"beta":2.71986,"ct0":1.28161335,"ct1":0.011719438,"ct2":0.0000892639,"k":0.00011249900000000001,"maximumFrequency":800000,"minimumFrequency":500000},{"alpha":2.61053,"beta":2.49772,"ct0":1.010843873,"ct1":0.006141983,"ct2":0.0000611871,"k":2.23928e-7,"maximumFrequency":1200000,"minimumFrequency":800000}],"referenceVolumetricLosses":null}]}},"numberStacks":1,"shape":{"aliases":["ELP 18/4/10","E 18/4/10/R","E18/4","E18/8"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01835,"minimum":0.017650000000000002,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0041,"minimum":0.0039000000000000003,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0102,"minimum":0.0098,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0021000000000000003,"minimum":0.0019,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0143,"minimum":0.0137,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0041,"minimum":0.0039000000000000003,"nominal":null}},"family":"planar e","familySubtype":null,"magneticCircuit":"open","name":"E 18/4/10/R","type":"standard"},"type":"two-piece set"},"geometricalDescription":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Ferroxcube","reference":"E18/4/10/R-3F36","status":"production"},"name":"E 18/4/10/R - 3F36 - Ungapped","processedDescription":{"columns":[{"area":0.00004,"coordinates":[0,0,0],"depth":0.01,"height":0.004,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.004},{"area":0.000021,"coordinates":[0.008,0,0],"depth":0.01,"height":0.004,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"lateral","width":0.002001},{"area":0.000021,"coordinates":[-0.008,0,0],"depth":0.01,"height":0.004,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"lateral","width":0.002001}],"depth":0.01,"effectiveParameters":{"effectiveArea":0.00004,"effectiveLength":0.024283185307179586,"effectiveVolume":9.713274122871836e-7,"minimumArea":0.00004},"height":0.008,"width":0.018000000000000002,"windingWindows":[{"angle":null,"area":0.00002,"coordinates":[0.002,0],"height":0.004,"radialHeight":null,"width":0.005}]}})");

        auto magneticFluxDensity = core.get_magnetic_flux_density_saturation(25.0, false);
        double expectedMagneticFluxDensity = 0.52;

        CHECK_CLOSE(expectedMagneticFluxDensity, magneticFluxDensity, expectedMagneticFluxDensity * maximumError);

        magneticFluxDensity = core.get_magnetic_flux_density_saturation(125.0, false);
        expectedMagneticFluxDensity = 0.4;
        CHECK_CLOSE(expectedMagneticFluxDensity, magneticFluxDensity, expectedMagneticFluxDensity * maximumError);
    }

    TEST(Test_Core_Functional_Description_Web_10) {
        auto coreJson = json::parse(R"({"distributorsInfo":[],"functionalDescription":{"coating":null,"gapping":[{"area":0.00003,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.007198,"distanceClosestParallelSurface":0.00435,"length":0.000005,"sectionDimensions":[0.0057,0.0051],"shape":"rectangular","type":"residual"},{"area":0.000015,"coordinates":[0.008625,0,0],"distanceClosestNormalSurface":0.007198,"distanceClosestParallelSurface":0.00435,"length":0.000005,"sectionDimensions":[0.00285,0.0051],"shape":"rectangular","type":"residual"},{"area":0.000015,"coordinates":[-0.008625,0,0],"distanceClosestNormalSurface":0.007198,"distanceClosestParallelSurface":0.00435,"length":0.000005,"sectionDimensions":[0.00285,0.0051],"shape":"rectangular","type":"residual"}],"material":"3C90","numberStacks":1,"shape":{"aliases":["R 140/103/25"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.14},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.103},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.025}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 140/103/25","type":"standard"},"type":"two-piece set"},"geometricalDescription":null,"manufacturerInfo":{"cost":null,"datasheetUrl":"https://ferroxcube.com/upload/media/product/file/Pr_ds/E20_10_5.pdf","family":null,"name":"Ferroxcube","orderCode":null,"reference":"E20/10/5-3C90","status":"production"},"name":"E 20/10/5 - 3C90 - Ungapped","processedDescription":null})");
        Core core(coreJson, true);
    }

    TEST(Missing_Core_Hermes) {
        // Check for segmentation fault
        auto coreJson = json::parse(R"({"functionalDescription": {"gapping": [], "material": "3C91", "numberStacks": 1, "shape": "P 11/7/I", "type": "two-piece set"}, "name": "temp"})");
        Core core(coreJson, true);

        auto functionalDescription = core.get_functional_description();
    }

    TEST(Test_Core_Initial_Permeability) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);

        Core core(coreJson, true);

        auto initialPermeability = core.get_initial_permeability(25);
        double expectedInitialPermeability = 2270;

        CHECK_CLOSE(expectedInitialPermeability, initialPermeability, expectedInitialPermeability * maximumError);

        initialPermeability = core.get_initial_permeability(125);
        expectedInitialPermeability = 3975;
        CHECK_CLOSE(expectedInitialPermeability, initialPermeability, expectedInitialPermeability * maximumError);
    }

    TEST(Test_Core_Effective_Permeability) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);

        Core core(coreJson, true);

        auto effectivePermeability = core.get_effective_permeability(25);
        double expectedEffectivePermeability = 136;

        CHECK_CLOSE(expectedEffectivePermeability, effectivePermeability, expectedEffectivePermeability * maximumError);

        effectivePermeability = core.get_effective_permeability(125);
        expectedEffectivePermeability = 139;
        CHECK_CLOSE(expectedEffectivePermeability, effectivePermeability, expectedEffectivePermeability * maximumError);
    }

    TEST(Test_Core_Reluctance) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);

        Core core(coreJson, true);

        auto reluctance = core.get_reluctance(25);
        double expectedReluctance = 1.02e+06;

        CHECK_CLOSE(expectedReluctance, reluctance, expectedReluctance * maximumError);

        reluctance = core.get_reluctance(125);
        expectedReluctance = 997019;
        CHECK_CLOSE(expectedReluctance, reluctance, expectedReluctance * maximumError);
    }

    TEST(Test_Core_Resistivity) {
        auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
        std::ifstream json_file(coreFilePath);

        auto coreJson = json::parse(json_file);

        Core core(coreJson, true);

        auto resistivity = core.get_resistivity(25);
        double expectedResistivity = 8;

        CHECK_CLOSE(expectedResistivity, resistivity, expectedResistivity * maximumError);

        resistivity = core.get_resistivity(125);
        expectedResistivity = 8;
        CHECK_CLOSE(expectedResistivity, resistivity, expectedResistivity * maximumError);
    }
}