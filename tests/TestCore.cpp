#include <source_location>
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "support/Painter.h"
#include "support/Settings.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/ReluctanceNetwork.h"
#include "physical_models/Inductance.h"
#include "processors/Inputs.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace { 

auto masPath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("MAS/").string();
std::string filePath = std::source_location::current().file_name();
double maximumError = 0.05;

TEST_CASE("E_55_21", "[constructive-model][core][processed-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    Core core(coreJson, true);

    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_55_21_N97_additive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000353 * numberStacks, 0.000353 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.124, 0.124 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(4.4e-05 * numberStacks, 4.4e-05 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.00035 * numberStacks, 0.00035 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.037, 0.037 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.01015, 0.01015 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0172, 0.0172 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.021 * numberStacks, 0.021 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.00935, 0.00935 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.021 * numberStacks, 0.021 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.00935, 0.00935 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.021 * numberStacks, 0.021 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("E_55_28_21", "[constructive-model][core][processed-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_28_21_3C95_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_55_28_21_3C95_additive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000353 * numberStacks, 0.000353 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.124, 0.124 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(4.4e-05 * numberStacks, 4.4e-05 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.00035 * numberStacks, 0.00035 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.037, 0.037 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.01015, 0.01015 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0172, 0.0172 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.021 * numberStacks, 0.021 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.00935, 0.00935 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.021 * numberStacks, 0.021 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.00935, 0.00935 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.021 * numberStacks, 0.021 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("E_19_8_5", "[constructive-model][core][processed-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_19_8_5_N87_substractive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_19_8_5_N87_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000225 * numberStacks, 0.0000225 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0396, 0.0396 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000891 * numberStacks, 0.000000891 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000221 * numberStacks, 0.0000221 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0114, 0.0114 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.00475, 0.00475 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0048, 0.0048 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0048 * numberStacks, 0.0048 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.00235, 0.00235 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0048 * numberStacks, 0.0048 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.00235, 0.00235 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0048 * numberStacks, 0.0048 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("ETD_39_20_13", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_ETD_39_20_13_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ETD 39/20/13";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_ETD_39_20_13_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000125 * numberStacks, 0.000125 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0922, 0.0922 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000011500 * numberStacks, 0.000011500 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000123 * numberStacks, 0.000123 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0282, 0.0282 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.00825, 0.00825 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0128, 0.0128 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0128 * numberStacks, 0.0128 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0048, 0.0048 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0128 * numberStacks, 0.0128 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.0048, 0.0048 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0128 * numberStacks, 0.0128 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("ETD_19_14_8", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_ETD_19_14_8_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ETD 19/14/8";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_ETD_19_14_8_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000441 * numberStacks, 0.0000441 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0553, 0.0553 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000002440 * numberStacks, 0.000002440 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000395 * numberStacks, 0.0000395 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0184, 0.0184 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0034, 0.0034 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0076, 0.0076 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0076 * numberStacks, 0.0076 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.00255, 0.00255 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0076 * numberStacks, 0.0076 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.00255, 0.00255 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0076 * numberStacks, 0.0076 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("ETD_54_28_19", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_ETD_54_28_19_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ETD 54/28/19";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_ETD_54_28_19_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000280 * numberStacks, 0.000280 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.127, 0.127 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000035600 * numberStacks, 0.000035600 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000280 * numberStacks, 0.000280 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0396, 0.0396 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0104, 0.0104 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0193, 0.0193 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0193 * numberStacks, 0.0193 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0072, 0.0072 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0193 * numberStacks, 0.0193 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.0072, 0.0072 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0193 * numberStacks, 0.0193 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("ER_54_18_18", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_ER_54_18_18_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ER 54/18/18";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_ER_54_18_18_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000256 * numberStacks, 0.000256 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.090, 0.090 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000023000 * numberStacks, 0.000023000 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000252 * numberStacks, 0.000252 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0216, 0.0216 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.01025, 0.01025 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0183, 0.0183 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0183 * numberStacks, 0.0183 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.007, 0.007 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0183 * numberStacks, 0.0183 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.007, 0.007 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0183 * numberStacks, 0.0183 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("ER_18_3_10", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_ER_18_3_10_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ER 18/3/10";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_ER_18_3_10_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000302 * numberStacks, 0.0000302 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0221, 0.0221 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000667 * numberStacks, 0.000000667 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000301 * numberStacks, 0.0000301 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0031, 0.0031 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0047, 0.0047 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0062, 0.0062 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0062 * numberStacks, 0.0062 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_minimum_width().value(), Catch::Matchers::WithinAbs(0.0012, 0.0012 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.01 * numberStacks, 0.01 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_minimum_width().value(), Catch::Matchers::WithinAbs(0.0012, 0.0012 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.01 * numberStacks, 0.01 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("E_102_20_38", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_E_102_20_38_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ELP 102/20/38";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 2;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_102_20_38_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000538 * numberStacks, 0.000538 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.1476, 0.1476 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000079410 * numberStacks, 0.000079410 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0005245 * numberStacks, 0.0005245 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0266, 0.0266 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.036, 0.036 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.014, 0.014 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0375 * numberStacks, 0.0375 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.008, 0.008 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0375 * numberStacks, 0.0375 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.008, 0.008 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0375 * numberStacks, 0.0375 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("E_14_3_5_5", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_E_14_3.5_5_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "ELP 14/3.5/5";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 3;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_14_3.5_5_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000143 * numberStacks, 0.0000143 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0207, 0.0207 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000296 * numberStacks, 0.000000296 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000139 * numberStacks, 0.0000139 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.004, 0.004 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.004, 0.004 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.003, 0.003 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.005 * numberStacks, 0.005 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0015, 0.0015 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.005 * numberStacks, 0.005 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.0015, 0.0015 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.005 * numberStacks, 0.005 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("EL_25_4_3", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_E_25_4.3_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EL 25/4.3";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_25_4.3_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000856 * numberStacks, 0.0000856 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.030, 0.030 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000002570 * numberStacks, 0.000002570 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000083 * numberStacks, 0.000083 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.004, 0.004 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.007255, 0.007255 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.00632, 0.00632 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.01454 * numberStacks, 0.01454 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.002085, 0.002085 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.002085, 0.002085 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("EL_11_2", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_E_11_2_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EL 11/2";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_E_11_2_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000165 * numberStacks, 0.0000165 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0137, 0.0137 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000226 * numberStacks, 0.000000226 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000159 * numberStacks, 0.0000159 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.002, 0.002 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.003195, 0.003195 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.00278, 0.00278 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0064 * numberStacks, 0.0064 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.000915, 0.000915 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0088 * numberStacks, 0.0088 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.000915, 0.000915 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0088 * numberStacks, 0.0088 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("EC_70", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EC_70_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EC 70";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EC_70_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000280 * numberStacks, 0.000280 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.144, 0.144 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000040420 * numberStacks, 0.000040420 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000211 * numberStacks, 0.000211 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0455, 0.0455 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0141, 0.0141 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0164, 0.0164 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0164 * numberStacks, 0.0164 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.013, 0.013 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0164 * numberStacks, 0.0164 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.013, 0.013 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0164 * numberStacks, 0.0164 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("EFD_10_5_3", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EFD_10_5_3_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EFD 10/5/3";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EFD_10_5_3_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000072 * numberStacks, 0.0000072 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0231, 0.0231 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000166 * numberStacks, 0.000000166 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000065 * numberStacks, 0.0000065 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0075, 0.0075 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.00155, 0.00155 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.00455, 0.00455 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.00145 * numberStacks, 0.00145 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.001425, 0.001425 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0027 * numberStacks, 0.0027 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.001425, 0.001425 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0027 * numberStacks, 0.0027 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("EFD_30_15_9", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EFD_30_15_9_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EFD 30/15/9";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EFD_30_15_9_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000069 * numberStacks, 0.000069 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.068, 0.068 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000004690 * numberStacks, 0.000004690 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000069 * numberStacks, 0.000069 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0224, 0.0224 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0039, 0.0039 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0146, 0.0146 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0049 * numberStacks, 0.0049 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0038, 0.0038 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0091 * numberStacks, 0.0091 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.0038, 0.0038 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0091 * numberStacks, 0.0091 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("EQ_30_8_20", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EQ_30_8_20_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EQ 30/8/20";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EQ_30_8_20_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000108 * numberStacks, 0.000108 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.046, 0.046 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000004970 * numberStacks, 0.000004970 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000095 * numberStacks, 0.000095 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0106, 0.0106 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0075, 0.0075 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.011, 0.011 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.011 * numberStacks, 0.011 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_minimum_width().value(), Catch::Matchers::WithinAbs(0.002, 0.002 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_minimum_width().value(), Catch::Matchers::WithinAbs(0.002, 0.002 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("EPX_10", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EPX_10_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EPX 10";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EPX_10_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000159 * numberStacks, 0.0000159 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0217, 0.0217 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000345 * numberStacks, 0.000000345 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000132 * numberStacks, 0.0000132 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0072, 0.0072 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.002825, 0.002825 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 3u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.00345, 0.00345 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.00485 * numberStacks, 0.00485 * numberStacks * 0.2));
    // Lateral depth updated after fixing the integer-division bug that dropped
    // the round bore term from the EP-family lateral leg area (CorePiece.cpp)
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.01027 * numberStacks, 0.01027 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.01027 * numberStacks, 0.01027 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("EPX_7", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EPX_7_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EPX 7";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EPX_7_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000172 * numberStacks, 0.0000172 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0157, 0.0157 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000270 * numberStacks, 0.000000270 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000139 * numberStacks, 0.0000139 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0045, 0.0045 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0019, 0.0019 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0034, 0.0034 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0057 * numberStacks, 0.0057 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0014 * numberStacks, 0.0014 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::OBLONG);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("EPO_13", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EPO_13_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EPO 13";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EPO_13_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000193 * numberStacks, 0.0000193 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0258, 0.0258 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000498 * numberStacks, 0.000000498 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000149 * numberStacks, 0.0000149 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.009, 0.009 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0026, 0.0026 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 3u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0045, 0.0045 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0045 * numberStacks, 0.0045 * numberStacks * 0.2));
    // Lateral depth updated after fixing the integer-division bug that dropped
    // the round bore term from the EP-family lateral leg area (CorePiece.cpp)
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.01068 * numberStacks, 0.01068 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.01068 * numberStacks, 0.01068 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("LP_42_25_15_8", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_LP_42_25_15_8_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "LP 42/25/15.8";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_LP_42_25_15_8_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000206 * numberStacks, 0.000206 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0901, 0.0901 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000018560 * numberStacks, 0.000018560 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000206 * numberStacks, 0.000206 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0228, 0.0228 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0095, 0.0095 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 3u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0162, 0.0162 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0162 * numberStacks, 0.0162 * numberStacks * 0.2));
    // Lateral depth updated after fixing the integer-division bug that dropped
    // the round bore term from the EP-family lateral leg area (CorePiece.cpp)
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0359 * numberStacks, 0.0359 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.0359 * numberStacks, 0.0359 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("EP_7", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EP_7_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EP 7";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EP_7_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000103 * numberStacks, 0.0000103 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0157, 0.0157 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000162 * numberStacks, 0.000000162 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000085 * numberStacks, 0.0000085 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.005, 0.005 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0019, 0.0019 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0034, 0.0034 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0034 * numberStacks, 0.0034 * numberStacks * 0.2));
    // Lateral width updated after fixing the integer-division bug that dropped
    // the round bore term from the EP-family lateral leg area (CorePiece.cpp)
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0256, 0.0256 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.00095 * numberStacks, 0.00095 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("EP_20", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_EP_20_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "EP 20";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_EP_20_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000078 * numberStacks, 0.000078 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.040, 0.040 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000003120 * numberStacks, 0.000003120 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000060 * numberStacks, 0.000060 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.014, 0.014 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.00355, 0.00355 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.009, 0.009 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.009 * numberStacks, 0.009 * numberStacks * 0.2));
    // Lateral width updated after fixing the integer-division bug that dropped
    // the round bore term from the EP-family lateral leg area (CorePiece.cpp)
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0807, 0.0807 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.00275 * numberStacks, 0.00275 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("RM_14", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_RM_14_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "RM 14";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_RM_14_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000200 * numberStacks, 0.000200 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.070, 0.070 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000014000 * numberStacks, 0.000014000 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000170 * numberStacks, 0.000170 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0208, 0.0208 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.007, 0.007 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.015, 0.015 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.015 * numberStacks, 0.015 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0066 * numberStacks, 0.0066 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.0066 * numberStacks, 0.0066 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("RM_7LP", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_RM_7LP_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "RM 7LP";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_RM_7LP_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    // ABT #783: these four were harvested in 2022, when "RM 7LP" was only an ALIAS pointing at
    // RM 7/10 — the WITH-centre-hole part. They were therefore RM 7's figures (Ae 40 mm2,
    // Ve 1190 mm3, Amin 32.3 mm2), and Amin went red the moment MAS 92a2af8 dropped the wrong LP
    // aliases and gave RM 7LP its own record. TDK rm_7.pdf, ordering code B65819P, header
    // "Without center hole", publishes for the LP set:
    //     le = 23.5 mm, Ae = 45.3 mm2, Amin = 39.6 mm2, Ve = 1060 mm3
    // MKF computes le 24.40 (+3.8 %), Ae 43.70 (-3.5 %), Amin 39.592 (-0.02 %), Ve 1066 (+0.6 %).
    // H is absent from the record BY DESIGN (no centre hole), so d4 = 0 is correct, not a
    // fallback; the subtype-2 branch of CorePieceRm never reads C, so its absence is inert.
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000453 * numberStacks, 0.0000453 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0235, 0.0235 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000001060 * numberStacks, 0.000001060 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000396 * numberStacks, 0.0000396 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0047, 0.0047 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.00375, 0.00375 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.00725, 0.00725 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.00725 * numberStacks, 0.00725 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.002775 * numberStacks, 0.002775 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.002775 * numberStacks, 0.002775 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("PQ_20_16", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_PQ_20_16_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "PQ 20/16";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_PQ_20_16_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000632 * numberStacks, 0.0000632 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0372, 0.0372 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000002360 * numberStacks, 0.000002360 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000544 * numberStacks, 0.0000544 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.01030, 0.01030 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0046, 0.0046 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0088, 0.0088 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0088 * numberStacks, 0.0088 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.014 * numberStacks, 0.014 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.014 * numberStacks, 0.014 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("PQ_107_87", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_PQ_107_87_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "PQ 107/87";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_PQ_107_87_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.001428 * numberStacks, 0.001428 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.204, 0.204 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000290600 * numberStacks, 0.000290600 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.001320 * numberStacks, 0.001320 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.056, 0.056 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.02635, 0.02635 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.041, 0.041 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.041 * numberStacks, 0.041 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.070 * numberStacks, 0.070 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_depth(), Catch::Matchers::WithinAbs(0.070 * numberStacks, 0.070 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("PM_114_93", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_PM_114_93_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "PM 114/93";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_PM_114_93_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.001720 * numberStacks, 0.001720 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.200, 0.200 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000344000 * numberStacks, 0.000344000 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.001380 * numberStacks, 0.001380 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.063, 0.063 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0225, 0.0225 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.043, 0.043 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.043 * numberStacks, 0.043 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.013 * numberStacks, 0.013 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.013 * numberStacks, 0.013 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("P_150_30", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_P_150_30_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "P 150/30";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_P_150_30_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.003580 * numberStacks, 0.003580 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.160, 0.160 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.00056600 * numberStacks, 0.00056600 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.002800 * numberStacks, 0.002800 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.030, 0.030 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0325, 0.0325 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.065, 0.065 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.065 * numberStacks, 0.065 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.010 * numberStacks, 0.010 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.010 * numberStacks, 0.00515 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("P_11_7", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_P_11_7_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "P 11/7";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_P_11_7_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0000162 * numberStacks, 0.0000162 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.0155, 0.0155 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000251 * numberStacks, 0.000000251 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.0000132 * numberStacks, 0.0000132 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.00440, 0.00440 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.00215, 0.00215 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0047, 0.0047 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0047 * numberStacks, 0.0047 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.00105 * numberStacks, 0.00105 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.00105 * numberStacks, 0.00105 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("P_7_4", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_P_7_4_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "P 7/4";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_P_7_4_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000007 * numberStacks, 0.000007 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.010, 0.010 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000070 * numberStacks, 0.000000070 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000006 * numberStacks, 0.000006 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0028, 0.0028 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0014, 0.0014 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.003, 0.003 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.003 * numberStacks, 0.003 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.00075 * numberStacks, 0.00075 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[2].get_width(), Catch::Matchers::WithinAbs(0.00075 * numberStacks, 0.00075 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::IRREGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[2].get_shape() == ColumnShape::IRREGULAR);
}

TEST_CASE("U_79_129_31", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_U_79_129_31_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "U 79/129/31";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_U_79_129_31_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000693 * numberStacks, 0.000693 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.309, 0.309 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000214220 * numberStacks, 0.000214220 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000693 * numberStacks, 0.000693 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.085, 0.085 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.034, 0.034 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.022, 0.022 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0315 * numberStacks, 0.0315 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.022 * numberStacks, 0.022 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0315 * numberStacks, 0.0315 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("U_26_22_16", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_U_26_22_16_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "U 26/22/16";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_U_26_22_16_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000131 * numberStacks, 0.000131 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.098, 0.098 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000012800 * numberStacks, 0.000012800 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000129 * numberStacks, 0.000129 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.026, 0.026 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.009, 0.009 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0084, 0.0084 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.016 * numberStacks, 0.016 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0084 * numberStacks, 0.0084 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.016 * numberStacks, 0.016 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("UR_48_39_17", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_UR_48_39_17_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "UR 48/39/17";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_UR_48_39_17_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000215 * numberStacks, 0.000215 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.186, 0.186 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000039990 * numberStacks, 0.000039990 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000215 * numberStacks, 0.000215 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0538, 0.0538 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0174, 0.0174 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.017, 0.017 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.017 * numberStacks, 0.017 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.013 * numberStacks, 0.013 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.017 * numberStacks, 0.017 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("UR_70_33_17", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_UR_70_33_17_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "UR 70/33/17";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_UR_70_33_17_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000214 * numberStacks, 0.000214 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.197, 0.197 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000043800 * numberStacks, 0.000043800 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000214 * numberStacks, 0.000214 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.0381, 0.0381 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.035, 0.035 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.01725, 0.01725 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.01725 * numberStacks, 0.01725 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.01725 * numberStacks, 0.01725 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.01725 * numberStacks, 0.01725 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::ROUND);
}

TEST_CASE("UR_55_39_36", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_UR_55_39_36_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "UR 55/38/36";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_UR_55_39_36_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000418 * numberStacks, 0.000418 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.188, 0.188 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000078570 * numberStacks, 0.000078570 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000418 * numberStacks, 0.000418 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.051, 0.051 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0196, 0.0196 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0235, 0.0235 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0235 * numberStacks, 0.0235 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.012 * numberStacks, 0.012 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.036 * numberStacks, 0.036 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("UR_64_40_20", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_UR_64_40_20_N97_substractive";
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "UR 64/40/20";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    REQUIRE(*(core.get_name()) == "core_UR_64_40_20_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000290 * numberStacks, 0.000290 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.210, 0.210 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000061000 * numberStacks, 0.000061000 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000290 * numberStacks, 0.000290 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.053, 0.053 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0232, 0.0232 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.020, 0.020 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.020 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::ROUND);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() == ColumnShape::ROUND);
}

TEST_CASE("UT_20", "[constructive-model][core][processed-description][smoke-test]") {
    json coreJson;
    coreJson["functionalDescription"] = json();
    coreJson["name"] = "core_UT_20_N97";
    coreJson["functionalDescription"]["type"] = "closedShape";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "UT 20";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson, true);
    double numberStacks = coreJson["functionalDescription"]["numberStacks"];

    auto geometrical_description = *(core.get_geometrical_description());

    REQUIRE(*(core.get_name()) == "core_UT_20_N97");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::CLOSED);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000013 * numberStacks, 0.000013 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.053, 0.053 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000000688 * numberStacks, 0.000000688 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000013 * numberStacks, 0.000013 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_height()), Catch::Matchers::WithinAbs(0.016, 0.016 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_width()), Catch::Matchers::WithinAbs(0.0075, 0.0075 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.0041, 0.0041 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.0046 * numberStacks, 0.0046 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_width(), Catch::Matchers::WithinAbs(0.0033 * numberStacks, 0.0033 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[1].get_depth(), Catch::Matchers::WithinAbs(0.0046 * numberStacks, 0.0046 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() ==
          ColumnShape::RECTANGULAR);
    REQUIRE(core.get_processed_description()->get_columns()[1].get_shape() ==
          ColumnShape::RECTANGULAR);
}

TEST_CASE("T_40_24_16", "[constructive-model][core][processed-description][smoke-test]") {
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

    REQUIRE(*(core.get_name()) == "core_T_40_24_16_N97_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.000125 * numberStacks, 0.000125 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_length(), Catch::Matchers::WithinAbs(0.09629, 0.09629 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_volume(), Catch::Matchers::WithinAbs(0.000012070 * numberStacks, 0.000012070 * numberStacks * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_minimum_area(), Catch::Matchers::WithinAbs(0.000125 * numberStacks, 0.000125 * numberStacks * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_radial_height()), Catch::Matchers::WithinAbs(0.012, 0.012 * 0.2));
    REQUIRE_THAT(*(core.get_processed_description()->get_winding_windows()[0].get_angle()), Catch::Matchers::WithinAbs(360, 360 * 0.2));
    REQUIRE(core.get_processed_description()->get_columns().size() == 1u);
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_width(), Catch::Matchers::WithinAbs(0.008, 0.008 * 0.2));
    REQUIRE_THAT(core.get_processed_description()->get_columns()[0].get_depth(), Catch::Matchers::WithinAbs(0.016 * numberStacks, 0.020 * numberStacks * 0.2));
    REQUIRE(core.get_processed_description()->get_columns()[0].get_shape() == ColumnShape::RECTANGULAR);
}

TEST_CASE("Core_Processed_Description_Web_0", "[constructive-model][core][processed-description][bug][smoke-test]") {
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

    REQUIRE(core.get_processed_description()->get_columns().size() == 2u);
}

TEST_CASE("Core_Processed_Description_Web_1", "[constructive-model][core][processed-description][bug][smoke-test]") {
    auto coreJson = json::parse(
        "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [], \"material\": \"3C97\", "
        "\"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": {\"A\": 0.0308, \"B\": 0.0264, \"C\": "
        "0.0265, \"D\": 0.016, \"E\": 0.01, \"G\": 0.0, \"H\": 0.0}, \"family\": \"u\", \"familySubtype\": \"1\", "
        "\"name\": \"Custom\", \"type\": \"custom\"}, \"type\": \"two-piece set\"}}");

    Core core(coreJson, true);

    REQUIRE_THAT(core.get_processed_description()->get_effective_parameters().get_effective_area(), Catch::Matchers::WithinAbs(0.0002756, 0.0002756 * 0.2));
    auto functionalDescription = core.get_functional_description();
}

TEST_CASE("Test_Core_All_Shapes", "[constructive-model][core][processed-description][smoke-test]") {
    settings.set_use_toroidal_cores(true);
    auto shapeNames = get_core_shape_names();
    for (auto shapeName : shapeNames) {
        if (shapeName.contains("PQI") || shapeName.contains("UI ")) {
            continue;
        }
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
        // if ((core.get_processed_description()->get_effective_parameters().get_effective_area() <= 0) ||
        //     (core.get_processed_description()->get_effective_parameters().get_effective_length() <= 0) ||
        //     (core.get_processed_description()->get_effective_parameters().get_effective_volume() <= 0) ||
        //     (core.get_processed_description()->get_effective_parameters().get_minimum_area() <= 0)) {

        //     std::cout << "shapeName: " << shapeName << std::endl;
        // }

        REQUIRE(core.get_processed_description()->get_effective_parameters().get_effective_area() > 0);
        REQUIRE(core.get_processed_description()->get_effective_parameters().get_effective_length() > 0);
        REQUIRE(core.get_processed_description()->get_effective_parameters().get_effective_volume() > 0);
        REQUIRE(core.get_processed_description()->get_effective_parameters().get_minimum_area() > 0);
    }
}

TEST_CASE("Core_Processed_Description_Web_2", "[constructive-model][core][processed-description][bug][smoke-test]") {
    auto coreJson = json::parse(R"({"functionalDescription": {"type": "twoPieceSet", "name": "150-2646", "shape": {"aliases": ["ER 9.5/5", "ER 9.5"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00955, "minimum": 0.00915, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00253, "minimum": 0.00238, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.005, "minimum": 0.0048, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00175, "minimum": 0.0016, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00775, "minimum": 0.0075, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00355, "minimum": 0.00325, "nominal": null}, "G": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0074, "minimum": 0.007, "nominal": null}}, "family": "planarER", "familySubtype": null, "magneticCircuit": "open", "name": "ER 9.5/2.5/5", "type": "standard"}, "material": {"alternatives": null, "application": null, "bhCycle": null, "coerciveForce": [{"magneticField": 27.2, "magneticFluxDensity": 0.0, "temperature": 100.0}, {"magneticField": 36.5, "magneticFluxDensity": 0.0, "temperature": 25.0}], "commercialName": null, "curieTemperature": 240.0, "density": 4700.0, "family": "TP", "heatCapacity": null, "heatConductivity": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "description": null, "family": null, "name": "TDG", "orderCode": null, "reference": null, "status": null}, "massLosses": null, "material": "ferrite", "materialComposition": "MnZn", "name": "TP5", "permeability": {"amplitude": null, "complex": null, "initial": [{"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -38.0172, "tolerance": null, "value": 964.291}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -32.5809, "tolerance": null, "value": 976.475}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -27.3041, "tolerance": null, "value": 989.988}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -22.0269, "tolerance": null, "value": 1007.57}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -16.7495, "tolerance": null, "value": 1026.22}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -11.6628, "tolerance": null, "value": 1053.38}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -6.19159, "tolerance": null, "value": 1087.54}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": -0.741648, "tolerance": null, "value": 1116.03}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 4.36765, "tolerance": null, "value": 1158.76}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 9.64787, "tolerance": null, "value": 1198.91}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 14.9285, "tolerance": null, "value": 1242.43}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 20.2095, "tolerance": null, "value": 1288.29}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 25.4907, "tolerance": null, "value": 1335.71}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 30.772100000000002, "tolerance": null, "value": 1384.68}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 36.3101, "tolerance": null, "value": 1440.25}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 42.6696, "tolerance": null, "value": 1503.07}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 47.6867, "tolerance": null, "value": 1537.3}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 52.8575, "tolerance": null, "value": 1586.52}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 58.1382, "tolerance": null, "value": 1630.56}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 63.4186, "tolerance": null, "value": 1671.75}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 68.6985, "tolerance": null, "value": 1709.57}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 73.9781, "tolerance": null, "value": 1744.79}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 79.3364, "tolerance": null, "value": 1771.53}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 84.5358, "tolerance": null, "value": 1804.85}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 89.8141, "tolerance": null, "value": 1830.21}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 95.0919, "tolerance": null, "value": 1852.2}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 100.369, "tolerance": null, "value": 1864.85}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 105.644, "tolerance": null, "value": 1866.34}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 110.918, "tolerance": null, "value": 1865.76}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 116.193, "tolerance": null, "value": 1862.84}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 121.387, "tolerance": null, "value": 1858.06}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 126.742, "tolerance": null, "value": 1855.7}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 132.016, "tolerance": null, "value": 1850.71}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 137.291, "tolerance": null, "value": 1848.05}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 142.565, "tolerance": null, "value": 1844.09}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 147.84, "tolerance": null, "value": 1842.99}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 153.115, "tolerance": null, "value": 1842.93}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 158.55, "tolerance": null, "value": 1843.81}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 163.665, "tolerance": null, "value": 1849.28}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 168.941, "tolerance": null, "value": 1855.18}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 174.217, "tolerance": null, "value": 1861.35}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 179.493, "tolerance": null, "value": 1868.03}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 184.769, "tolerance": null, "value": 1875.23}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 190.044, "tolerance": null, "value": 1880.61}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 195.319, "tolerance": null, "value": 1880.29}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 200.593, "tolerance": null, "value": 1874.25}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 207.786, "tolerance": null, "value": 1871.05}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 213.062, "tolerance": null, "value": 1879.29}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 218.338, "tolerance": null, "value": 1889.6}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 223.614, "tolerance": null, "value": 1898.88}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 228.89, "tolerance": null, "value": 1906.34}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 234.166, "tolerance": null, "value": 1914.06}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 239.521, "tolerance": null, "value": 1917.97}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 244.719, "tolerance": null, "value": 1933.65}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 249.996, "tolerance": null, "value": 1952.0}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 255.274, "tolerance": null, "value": 1978.06}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 260.707, "tolerance": null, "value": 1965.15}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 263.777, "tolerance": null, "value": 1911.43}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 265.92, "tolerance": null, "value": 1491.87}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 264.844, "tolerance": null, "value": 1819.13}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 266.953, "tolerance": null, "value": 994.433}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 265.723, "tolerance": null, "value": 1733.09}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 265.973, "tolerance": null, "value": 1652.04}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 266.251, "tolerance": null, "value": 1577.83}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 266.707, "tolerance": null, "value": 1402.02}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 266.539, "tolerance": null, "value": 1337.81}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 267.166, "tolerance": null, "value": 1247.08}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 266.996, "tolerance": null, "value": 1170.84}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 267.338, "tolerance": null, "value": 1094.92}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 268.075, "tolerance": null, "value": 867.595}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 267.904, "tolerance": null, "value": 781.74}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 268.452, "tolerance": null, "value": 698.017}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 268.363, "tolerance": null, "value": 624.291}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 268.993, "tolerance": null, "value": 554.361}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 269.313, "tolerance": null, "value": 383.526}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 269.278, "tolerance": null, "value": 295.597}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 269.158, "tolerance": null, "value": 471.599}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 269.798, "tolerance": null, "value": 193.482}, {"frequency": null, "magneticFieldDcBias": null, "magneticFluxDensityPeak": null, "modifiers": null, "temperature": 269.896, "tolerance": null, "value": 130.059}]}, "remanence": [{"magneticField": 0.0, "magneticFluxDensity": 0.098, "temperature": 100.0}, {"magneticField": 0.0, "magneticFluxDensity": 0.14, "temperature": 25.0}], "resistivity": [{"temperature": 20.0, "value": 8.0}], "saturation": [{"magneticField": 1194.0, "magneticFluxDensity": 0.38, "temperature": 100.0}, {"magneticField": 1194.0, "magneticFluxDensity": 0.47000000000000003, "temperature": 25.0}], "type": "commercial", "volumetricLosses": {"default": [{"a": null, "b": null, "c": null, "coefficients": null, "d": null, "factors": null, "method": "steinmetz", "ranges": [{"alpha": 1.5832445907356367, "beta": 3.447276210669818, "ct0": 0.0006981954642389527, "ct1": -0.02648487944707231, "ct2": -0.0001649186139824856, "k": 2.0736187833042155, "maximumFrequency": 150000.0, "minimumFrequency": 1.0}, {"alpha": 2.0382319724990148, "beta": 3.044361726422575, "ct0": 10.063670862549936, "ct1": 7.787255730814788, "ct2": -101.06805047707402, "k": 0.002087675611586363, "maximumFrequency": 1000000.0, "minimumFrequency": 150000.0}, {"alpha": 1.443736158982313, "beta": 2.492266808268583, "ct0": 0.0006461213260757808, "ct1": -0.024170654037881944, "ct2": -0.00014177115251103593, "k": 1.7406282450344872, "maximumFrequency": 1000000000.0, "minimumFrequency": 1000000.0}], "referenceVolumetricLosses": null}]}}, "numberStacks": 1, "gapping": []}})");

    Core ungappedCore(coreJson, true);

    OpenMagnetics::Coil coil;
    OpenMagnetics::Winding winding;
    winding.set_number_turns(1);
    coil.get_mutable_functional_description().push_back(winding);

    OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
    ungappedCore.process_data();
    OpenMagnetics::Inputs inputs;
    MAS::DesignRequirements designRequirements;
    DimensionWithTolerance inductanceWithTolerance;
    inductanceWithTolerance.set_nominal(OpenMagnetics::roundFloat(4.5000000000000003e-07, 10));
    designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
    inputs.set_design_requirements(designRequirements);

    auto gapping =  magnetizingInductanceModel.calculate_gapping_from_number_turns_and_inductance(ungappedCore, coil, &inputs, OpenMagnetics::GappingType::GROUND);
    // Repro point: a ground gap must be computable for the requested inductance.
    REQUIRE(gapping.size() > 0);
    for (auto& gap : gapping) {
        CHECK(gap.get_length() > 0);
    }
    ungappedCore.get_mutable_functional_description().set_gapping(gapping);
    ungappedCore.process_gap();
    auto functionalDescription = ungappedCore.get_functional_description();
    REQUIRE(functionalDescription.get_gapping().size() == gapping.size());
}

TEST_CASE("E_19_8_5_Geometrical_Description", "[constructive-model][core][geometrical-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_19_8_5_N87_substractive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    Core core(coreJson, true);

    auto geometrical_description = *(core.get_geometrical_description());

    REQUIRE(*(core.get_name()) == "core_E_19_8_5_N87_substractive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE(geometrical_description.size() == 2u);
    REQUIRE(geometrical_description[0].get_machining());
    REQUIRE(!geometrical_description[1].get_machining());
    REQUIRE(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
    REQUIRE(geometrical_description[1].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
}

TEST_CASE("E_55_21_GeometricalDescription", "[constructive-model][core][geometrical-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    Core core(coreJson, true);

    auto geometrical_description = *(core.get_geometrical_description());

    REQUIRE(*(core.get_name()) == "core_E_55_21_N97_additive");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE(geometrical_description.size() == 6u);
    REQUIRE(!geometrical_description[0].get_machining());
    REQUIRE(!geometrical_description[1].get_machining());
    REQUIRE(!geometrical_description[2].get_machining());
    REQUIRE(!geometrical_description[3].get_machining());
    REQUIRE(!geometrical_description[4].get_machining());
    REQUIRE(!geometrical_description[5].get_machining());
    REQUIRE(geometrical_description[4].get_type() == CoreGeometricalDescriptionElementType::SPACER);
    REQUIRE(geometrical_description[5].get_type() == CoreGeometricalDescriptionElementType::SPACER);
}

TEST_CASE("T_40_24_16_2", "[constructive-model][core][geometrical-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_T_40_24_16_N97.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    Core core(coreJson, true);

    auto geometrical_description = *(core.get_geometrical_description());

    REQUIRE(*(core.get_name()) == "core_T_40_24_16_N97");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE(geometrical_description.size() == 1u);
}

TEST_CASE("Core_Web_0", "[constructive-model][core][geometrical-description][smoke-test]") {
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
        "familySubtype": "1", "name": "custom", "type": "custom"}, "type": "twoPieceSet"}})");

    Core core(coreJson, true);

    auto geometrical_description = *(core.get_geometrical_description());

    REQUIRE(*(core.get_name()) == "Custom_0");
    REQUIRE(std::get<CoreMaterial>(core.get_mutable_functional_description().get_mutable_material())
              .get_mutable_volumetric_losses()["default"]
              .size() > 0);
    REQUIRE(geometrical_description.size() == 2u);
    REQUIRE(geometrical_description[0].get_machining());
    REQUIRE(geometrical_description[1].get_machining());
    REQUIRE(geometrical_description[0].get_machining()->size() == 2);
    REQUIRE(geometrical_description[1].get_machining()->size() == 2);
    REQUIRE(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
    REQUIRE(geometrical_description[1].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
}

TEST_CASE("Test_Core_Geometrical_Description_Web_1", "[constructive-model][core][geometrical-description][smoke-test]") {
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

    REQUIRE(geometrical_description.size() == 2u);
    REQUIRE(geometrical_description[0].get_machining());
    REQUIRE(geometrical_description[1].get_machining());
    REQUIRE(geometrical_description[0].get_machining()->size() == 3);
    REQUIRE(geometrical_description[1].get_machining()->size() == 2);
    REQUIRE(geometrical_description[0].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
    REQUIRE(geometrical_description[1].get_type() == CoreGeometricalDescriptionElementType::HALF_SET);
}

TEST_CASE("E_55_21_all_gaps_residual", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    coreJson["functionalDescription"]["gapping"][0] = coreJson["functionalDescription"]["gapping"][1];

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();

    REQUIRE(functionalDescription.get_gapping().size() == 3u);
    REQUIRE(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
    REQUIRE(*(functionalDescription.get_gapping()[0].get_shape()) ==
          *(functionalDescription.get_gapping()[1].get_shape()));
    REQUIRE(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) ==
          *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
    REQUIRE(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[1].get_area()) * 2, 0.2));
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[0] ==
          -(*functionalDescription.get_gapping()[2].get_coordinates())[0]);
}

TEST_CASE("E_55_21_central_gap", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();

    REQUIRE(functionalDescription.get_gapping().size() == 3u);
    REQUIRE(functionalDescription.get_gapping()[0].get_type() != functionalDescription.get_gapping()[1].get_type());
    REQUIRE(*(functionalDescription.get_gapping()[0].get_shape()) ==
          *(functionalDescription.get_gapping()[1].get_shape()));
    REQUIRE(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) !=
          *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
    REQUIRE(functionalDescription.get_gapping()[0].get_length() != functionalDescription.get_gapping()[1].get_length());
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[1].get_area()) * 2, 0.2));
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] != 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[0] ==
          -(*functionalDescription.get_gapping()[2].get_coordinates())[0]);
}

TEST_CASE("E_55_21_gap_all_columns", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    coreJson["functionalDescription"]["gapping"][1] = coreJson["functionalDescription"]["gapping"][0];
    coreJson["functionalDescription"]["gapping"][2] = coreJson["functionalDescription"]["gapping"][0];

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();

    REQUIRE(functionalDescription.get_gapping().size() == 3u);
    REQUIRE(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
    REQUIRE(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[2].get_type());
    REQUIRE(*(functionalDescription.get_gapping()[0].get_shape()) ==
          *(functionalDescription.get_gapping()[1].get_shape()));
    REQUIRE(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) ==
          *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
    REQUIRE(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[1].get_area()) * 2, 0.2));
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[0] ==
          -(*functionalDescription.get_gapping()[2].get_coordinates())[0]);
}

TEST_CASE("E_55_21_central_distributed_gap_even", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    coreJson["functionalDescription"]["gapping"].push_back(coreJson["functionalDescription"]["gapping"][0]);

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();

    REQUIRE(functionalDescription.get_gapping().size() == 4u);
    REQUIRE(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
    REQUIRE(functionalDescription.get_gapping()[0].get_type() != functionalDescription.get_gapping()[2].get_type());
    REQUIRE(*(functionalDescription.get_gapping()[0].get_shape()) ==
          *(functionalDescription.get_gapping()[1].get_shape()));
    REQUIRE(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) ==
          *(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()));
    REQUIRE(*(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()) !=
          *(functionalDescription.get_gapping()[2].get_distance_closest_normal_surface()));
    REQUIRE(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
    REQUIRE(functionalDescription.get_gapping()[0].get_length() != functionalDescription.get_gapping()[2].get_length());
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[1].get_area()), 0.2));
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[2].get_area()) * 2, 0.2));
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] ==
          (*functionalDescription.get_gapping()[1].get_coordinates())[0]);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] ==
          -(*functionalDescription.get_gapping()[1].get_coordinates())[1]);
    REQUIRE((*functionalDescription.get_gapping()[2].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[2].get_coordinates())[0] ==
          -(*functionalDescription.get_gapping()[3].get_coordinates())[0]);
}

TEST_CASE("E_55_21_central_distributed_gap_odd", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);
    coreJson["functionalDescription"]["gapping"].push_back(coreJson["functionalDescription"]["gapping"][0]);
    coreJson["functionalDescription"]["gapping"].push_back(coreJson["functionalDescription"]["gapping"][0]);

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();

    REQUIRE(functionalDescription.get_gapping().size() == 5u);
    REQUIRE(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[1].get_type());
    REQUIRE(functionalDescription.get_gapping()[0].get_type() == functionalDescription.get_gapping()[2].get_type());
    REQUIRE(functionalDescription.get_gapping()[0].get_type() != functionalDescription.get_gapping()[3].get_type());
    REQUIRE(*(functionalDescription.get_gapping()[0].get_shape()) ==
          *(functionalDescription.get_gapping()[1].get_shape()));
    REQUIRE(*(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()) >
          *(functionalDescription.get_gapping()[0].get_distance_closest_normal_surface()));
    REQUIRE(*(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()) >
          *(functionalDescription.get_gapping()[2].get_distance_closest_normal_surface()));
    REQUIRE(*(functionalDescription.get_gapping()[1].get_distance_closest_normal_surface()) <
          *(functionalDescription.get_gapping()[3].get_distance_closest_normal_surface()));
    REQUIRE(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[1].get_length());
    REQUIRE(functionalDescription.get_gapping()[0].get_length() == functionalDescription.get_gapping()[2].get_length());
    REQUIRE(functionalDescription.get_gapping()[0].get_length() != functionalDescription.get_gapping()[3].get_length());
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[1].get_area()), 0.2));
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[2].get_area()), 0.2));
    REQUIRE_THAT(*(functionalDescription.get_gapping()[0].get_area()), Catch::Matchers::WithinAbs(*(functionalDescription.get_gapping()[3].get_area()) * 2, 0.2));
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] ==
          (*functionalDescription.get_gapping()[1].get_coordinates())[0]);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] ==
          (*functionalDescription.get_gapping()[2].get_coordinates())[0]);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] ==
          -(*functionalDescription.get_gapping()[2].get_coordinates())[1]);
    REQUIRE((*functionalDescription.get_gapping()[3].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[3].get_coordinates())[0] ==
          -(*functionalDescription.get_gapping()[4].get_coordinates())[0]);
}

TEST_CASE("Test_Core_Functional_Description_Web_0", "[constructive-model][core][functional-description][bug][smoke-test]") {
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
    // Repro point: construction must fully process the core (these were segfault/crash repros).
    REQUIRE(core.get_processed_description());
    auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
    CHECK(effectiveParameters.get_effective_area() > 0);
    CHECK(effectiveParameters.get_effective_length() > 0);
    CHECK(effectiveParameters.get_effective_volume() > 0);
}

TEST_CASE("Test_Core_Functional_Description_Web_1", "[constructive-model][core][functional-description][bug][smoke-test]") {
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

    REQUIRE(functionalDescription.get_gapping().size() == 2u);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
    // ABT #644: this used to require [1] == 0. A single ground gap is machined into ONE half --
    // grinding one piece by the whole gap is cheaper than grinding two by half each -- so it spans
    // 0..length and its centre sits at +length/2, here 0.0001/2. The old expectation came from the
    // distributed branch, which a bare one-gap list used to fall through to; it contradicted
    // TestCore.cpp's own E_55_21_central_gap (which requires [1] != 0) and
    // E_19_8_5_Geometrical_Description (which requires exactly one machined half-set).
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0.0001 / 2);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[2] == 0);

    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[2] != 0);
}

TEST_CASE("Test_Core_Functional_Description_Web_2", "[constructive-model][core][functional-description][bug][smoke-test]") {
    // Tests that a distributed but aligned gapping does not get recalculated
    auto coreJson = json::parse(
        "{\"name\": \"default\", \"functionalDescription\": {\"gapping\": [{\"area\": 1.5e-05, \"coordinates\": "
        "[0.0, 0.0, 0.0], \"distanceClosestNormalSurface\": 0.0041, \"length\": 0.001, \"sectionDimensions\": "
        "[0.0043, 0.0043], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 1.5e-05, \"coordinates\": "
        "[0.0, 0.001, 0.0], \"distanceClosestNormalSurface\": 0.0041, \"length\": 0.001, \"sectionDimensions\": "
        "[0.0043, 0.0043], \"shape\": \"round\", \"type\": \"subtractive\"}, {\"area\": 8.8e-05, \"coordinates\": "
        "[0.0, 0.0, -0.005751], \"distanceClosestNormalSurface\": 0.004598, \"length\": 5e-06, "
        "\"sectionDimensions\": [0.058628, 0.001501], \"shape\": \"irregular\", \"type\": \"residual\"}"
        "], \"material\": \"3C97\", \"numberStacks\": 1, \"shape\": {\"aliases\": [], \"dimensions\": "
        "{\"A\": 0.0125, \"B\": 0.0064, \"C\": 0.0088, \"D\": 0.0046, \"E\": 0.01, \"F\": 0.0043, \"G\": 0.000, "
        "\"H\": 0.0, \"K\": 0.0023}, \"family\": \"ep\", \"familySubtype\": \"1\", \"name\": \"Custom\", \"type\": "
        "\"custom\"}, \"type\": \"two-piece set\"}, \"geometricalDescription\": null, \"processedDescription\": "
        "null}");

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();

    REQUIRE(functionalDescription.get_gapping().size() == 3u);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[0].get_coordinates())[2] == 0);

    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[1] != 0);
    REQUIRE((*functionalDescription.get_gapping()[1].get_coordinates())[2] == 0);

    REQUIRE((*functionalDescription.get_gapping()[2].get_coordinates())[0] == 0);
    REQUIRE((*functionalDescription.get_gapping()[2].get_coordinates())[1] == 0);
    REQUIRE((*functionalDescription.get_gapping()[2].get_coordinates())[2] != 0);
}

TEST_CASE("Test_Core_Functional_Description_Web_3", "[constructive-model][core][functional-description][bug][smoke-test]") {
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
    // Repro point: construction must fully process the core (these were segfault/crash repros).
    REQUIRE(core.get_processed_description());
    auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
    CHECK(effectiveParameters.get_effective_area() > 0);
    CHECK(effectiveParameters.get_effective_length() > 0);
    CHECK(effectiveParameters.get_effective_volume() > 0);
}

TEST_CASE("Test_Core_Functional_Description_Web_4", "[constructive-model][core][functional-description][bug][smoke-test]") {
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
    // Repro point: construction must fully process the core (these were segfault/crash repros).
    REQUIRE(core.get_processed_description());
    auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
    CHECK(effectiveParameters.get_effective_area() > 0);
    CHECK(effectiveParameters.get_effective_length() > 0);
    CHECK(effectiveParameters.get_effective_volume() > 0);
}

TEST_CASE("Test_Core_Functional_Description_Web_5", "[constructive-model][core][functional-description][bug][smoke-test]") {
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
    REQUIRE(functionalDescription.get_gapping().size() == 5u);
}

TEST_CASE("Test_Core_Functional_Description_Web_6", "[constructive-model][core][functional-description][bug][smoke-test]") {
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
    REQUIRE(functionalDescription.get_gapping().size() == 2u);
}

TEST_CASE("Test_Core_Functional_Description_Web_7", "[constructive-model][core][functional-description][bug][smoke-test]") {
    // Check for segmentation fault
    auto coreJson = json::parse(R"({"name":"My Core","functionalDescription":{"coating":null,"gapping":[{"area":0.000057,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.010097499999999999,"distanceClosestParallelSurface":0.005050000000000001,"length":0.000005,"sectionDimensions":[0.0085,0.0085],"shape":"round","type":"residual"},{"area":0.000028,"coordinates":[0.01075,0,0],"distanceClosestNormalSurface":0.010097499999999999,"distanceClosestParallelSurface":0.005050000000000001,"length":0.000005,"sectionDimensions":[0.0029,0.0085],"shape":"irregular","type":"residual"},{"area":0.000028,"coordinates":[-0.01075,0,0],"distanceClosestNormalSurface":0.010097499999999999,"distanceClosestParallelSurface":0.005050000000000001,"length":0.000005,"sectionDimensions":[0.0029,0.0085],"shape":"irregular","type":"residual"}],"material":"3C97","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":0.0576,"B":0.028399999999999998,"C":0.0155,"D":0.016,"H":0.0159,"G":0},"family":"ur","familySubtype":"2","name":"UR 57/28/16","type":"standard"},"type":"twoPieceSet"}})");

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();
    REQUIRE(functionalDescription.get_gapping().size() == 2u);
}

TEST_CASE("Test_Core_Functional_Description_Web_8", "[constructive-model][core][functional-description][bug][smoke-test]") {
    // Check for segmentation fault
    auto coreJson = json::parse(R"({"functionalDescription": {"type": "twoPieceSet", "material": "3C97", "shape": "U 80/150/30", "gapping": [{"length": 0.003, "type": "additive", "coordinates": [0, 0, 0 ] }, {"length": 0.003, "type": "additive", "coordinates": [0.0595, 0, 0 ] } ], "numberStacks": 1 }, "name": "My Core", "geometricalDescription": null, "processedDescription": null })");

    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();
    REQUIRE(functionalDescription.get_gapping().size() == 2u);
}

TEST_CASE("Test_Core_Functional_Description_Web_9", "[constructive-model][core][functional-description][bug][smoke-test]") {
    Core core = json::parse(R"({"distributorsInfo":[{"cost":1.17,"country":"USA","distributedArea":"International","email":null,"link":"https://www.digikey.com/en/products/detail/ferroxcube/E18-4-10-R-3F36/7041469","name":"Digi-Key","phone":null,"quantity":7063,"reference":"1779-1009-ND","updatedAt":"05/10/2023"},{"cost":null,"country":"UK","distributedArea":"International","email":null,"link":"https://www.shop.gatewaycando.com/magnetics/cores","name":"Gateway","phone":null,"quantity":5,"reference":"E18/4/10/R-3F36","updatedAt":"05/10/2023"}],"functionalDescription":{"coating":null,"gapping":[{"area":0.00004,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.0019975,"distanceClosestParallelSurface":0.005,"length":0.000005,"sectionDimensions":[0.004,0.01],"shape":"rectangular","type":"residual"},{"area":0.000021,"coordinates":[0.008,0,0],"distanceClosestNormalSurface":0.0019975,"distanceClosestParallelSurface":0.005,"length":0.000005,"sectionDimensions":[0.002001,0.01],"shape":"rectangular","type":"residual"},{"area":0.000021,"coordinates":[-0.008,0,0],"distanceClosestNormalSurface":0.0019975,"distanceClosestParallelSurface":0.005,"length":0.000005,"sectionDimensions":[0.002001,0.01],"shape":"rectangular","type":"residual"}],"material":{"bhCycle":null,"coerciveForce":[{"magneticField":32,"magneticFluxDensity":0,"temperature":100},{"magneticField":37,"magneticFluxDensity":0,"temperature":25}],"curieTemperature":230,"density":4750,"family":"3F","heatCapacity":{"excludeMaximum":null,"excludeMinimum":null,"maximum":800,"minimum":700,"nominal":null},"heatConductivity":{"excludeMaximum":null,"excludeMinimum":null,"maximum":5,"minimum":3.5,"nominal":null},"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Ferroxcube","reference":null,"status":null},"material":"ferrite","name":"3F36","permeability":{"amplitude":null,"initial":[{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-40,"tolerance":null,"value":1577},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-30,"tolerance":null,"value":1590},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-20,"tolerance":null,"value":1611},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":-10,"tolerance":null,"value":1633},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":0,"tolerance":null,"value":1657},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":10,"tolerance":null,"value":1683},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":20,"tolerance":null,"value":1710},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":30,"tolerance":null,"value":1767},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":40,"tolerance":null,"value":1792},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":50,"tolerance":null,"value":1818},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":60,"tolerance":null,"value":1842},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":70,"tolerance":null,"value":1780},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":80,"tolerance":null,"value":1794},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":90,"tolerance":null,"value":1807},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":100,"tolerance":null,"value":1818},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":110,"tolerance":null,"value":1893},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":120,"tolerance":null,"value":1897},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":130,"tolerance":null,"value":1901},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":140,"tolerance":null,"value":1907},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":150,"tolerance":null,"value":1848},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":160,"tolerance":null,"value":1855},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":170,"tolerance":null,"value":1865},{"frequency":10000,"magneticFieldDcBias":null,"magneticFluxDensityPeak":null,"modifiers":null,"temperature":180,"tolerance":null,"value":1878}]},"remanence":[{"magneticField":0,"magneticFluxDensity":0.105,"temperature":100},{"magneticField":0,"magneticFluxDensity":0.125,"temperature":25}],"resistivity":[{"temperature":25,"value":12}],"saturation":[{"magneticField":1200,"magneticFluxDensity":0.42,"temperature":100},{"magneticField":1200,"magneticFluxDensity":0.52,"temperature":25}],"type":"commercial","volumetricLosses":{"default":[{"a":null,"b":null,"c":null,"coefficients":{"excessLossesCoefficient":1.82280296e-20,"resistivityFrequencyCoefficient":4.25672064e-28,"resistivityMagneticFluxDensityCoefficient":14.7587264,"resistivityOffset":6.465016450000001e-17,"resistivityTemperatureCoefficient":5.84938089e-16},"d":null,"method":"roshen","ranges":null,"referenceVolumetricLosses":null},{"a":null,"b":null,"c":null,"coefficients":null,"d":null,"method":"steinmetz","ranges":[{"alpha":1.43902,"beta":3.26718,"ct0":1.232717265,"ct1":0.010783518,"ct2":0.00008394600000000001,"k":6.83,"maximumFrequency":499999,"minimumFrequency":100000},{"alpha":2.19515,"beta":2.71986,"ct0":1.28161335,"ct1":0.011719438,"ct2":0.0000892639,"k":0.00011249900000000001,"maximumFrequency":800000,"minimumFrequency":500000},{"alpha":2.61053,"beta":2.49772,"ct0":1.010843873,"ct1":0.006141983,"ct2":0.0000611871,"k":2.23928e-7,"maximumFrequency":1200000,"minimumFrequency":800000}],"referenceVolumetricLosses":null}]}},"numberStacks":1,"shape":{"aliases":["ELP 18/4/10","E 18/4/10/R","E18/4","E18/8"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01835,"minimum":0.017650000000000002,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0041,"minimum":0.0039000000000000003,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0102,"minimum":0.0098,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0021000000000000003,"minimum":0.0019,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0143,"minimum":0.0137,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0041,"minimum":0.0039000000000000003,"nominal":null}},"family":"planarE","familySubtype":null,"magneticCircuit":"open","name":"E 18/4/10/R","type":"standard"},"type":"twoPieceSet"},"geometricalDescription":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Ferroxcube","reference":"E18/4/10/R-3F36","status":"production"},"name":"E 18/4/10/R - 3F36 - Ungapped","processedDescription":{"columns":[{"area":0.00004,"coordinates":[0,0,0],"depth":0.01,"height":0.004,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.004},{"area":0.000021,"coordinates":[0.008,0,0],"depth":0.01,"height":0.004,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"lateral","width":0.002001},{"area":0.000021,"coordinates":[-0.008,0,0],"depth":0.01,"height":0.004,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"lateral","width":0.002001}],"depth":0.01,"effectiveParameters":{"effectiveArea":0.00004,"effectiveLength":0.024283185307179586,"effectiveVolume":9.713274122871836e-7,"minimumArea":0.00004},"height":0.008,"width":0.018000000000000002,"windingWindows":[{"angle":null,"area":0.00002,"coordinates":[0.002,0],"height":0.004,"radialHeight":null,"width":0.005}]}})");

    auto magneticFluxDensity = core.get_magnetic_flux_density_saturation(25.0, false);
    double expectedMagneticFluxDensity = 0.52;

    REQUIRE_THAT(expectedMagneticFluxDensity, Catch::Matchers::WithinAbs(magneticFluxDensity, expectedMagneticFluxDensity * maximumError));

    magneticFluxDensity = core.get_magnetic_flux_density_saturation(125.0, false);
    expectedMagneticFluxDensity = 0.4;
    REQUIRE_THAT(expectedMagneticFluxDensity, Catch::Matchers::WithinAbs(magneticFluxDensity, expectedMagneticFluxDensity * maximumError));
}

TEST_CASE("Test_Core_Functional_Description_Web_10", "[constructive-model][core][functional-description][bug]") {
    auto coreJson = json::parse(R"({"distributorsInfo":[],"functionalDescription":{"coating":null,"gapping":[{"area":0.00003,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.007198,"distanceClosestParallelSurface":0.00435,"length":0.000005,"sectionDimensions":[0.0057,0.0051],"shape":"rectangular","type":"residual"},{"area":0.000015,"coordinates":[0.008625,0,0],"distanceClosestNormalSurface":0.007198,"distanceClosestParallelSurface":0.00435,"length":0.000005,"sectionDimensions":[0.00285,0.0051],"shape":"rectangular","type":"residual"},{"area":0.000015,"coordinates":[-0.008625,0,0],"distanceClosestNormalSurface":0.007198,"distanceClosestParallelSurface":0.00435,"length":0.000005,"sectionDimensions":[0.00285,0.0051],"shape":"rectangular","type":"residual"}],"material":"3C90","numberStacks":1,"shape":{"aliases":["R 140/103/25"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.14},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.103},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.025}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 140/103/25","type":"standard"},"type":"twoPieceSet"},"geometricalDescription":null,"manufacturerInfo":{"cost":null,"datasheetUrl":"https://ferroxcube.com/upload/media/product/file/Pr_ds/E20_10_5.pdf","family":null,"name":"Ferroxcube","orderCode":null,"reference":"E20/10/5-3C90","status":"production"},"name":"E 20/10/5 - 3C90 - Ungapped","processedDescription":null})");
    SKIP("Test needs investigation");
    Core core(coreJson, true);
}

TEST_CASE("Missing_Core_Hermes", "[constructive-model][core][functional-description][bug][smoke-test]") {
    // Check for segmentation fault
    auto coreJson = json::parse(R"({"functionalDescription": {"gapping": [], "material": "3C91", "numberStacks": 1, "shape": "P 11/7/I", "type": "twoPieceSet"}, "name": "temp"})");
    Core core(coreJson, true);

    auto functionalDescription = core.get_functional_description();
    // Repro point: construction must fully process the core (these were segfault/crash repros).
    REQUIRE(core.get_processed_description());
    auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
    CHECK(effectiveParameters.get_effective_area() > 0);
    CHECK(effectiveParameters.get_effective_length() > 0);
    CHECK(effectiveParameters.get_effective_volume() > 0);
}

TEST_CASE("Test_Core_Standard_Shape_Reference_By_Object", "[constructive-model][core][functional-description][bug][smoke-test]") {
    // ABT #626: a core whose shape is a *reference* to the standard catalog given as a full
    // object ({"type": "standard", "family": <fam>, "name": <catalog name>}, no "dimensions" --
    // legal per core/shape.json, which does not require "dimensions") used to throw "bad
    // optional access" for EVERY family/name, 100% reproduction. Root cause: Core::process_data()
    // only resolved the database lookup when get_shape() held a bare std::string; when it held a
    // CoreShape object instead (this form), the lookup was skipped entirely and
    // CorePiece::factory() received a CoreShape with no "dimensions" map at all, so every
    // per-family builder's `get_shape().get_dimensions().value()` dereferenced a disengaged
    // optional. This mirrors the exact repro reported against PyOM's
    // calculate_core_processed_description (Core(json, false, false, false) +
    // core.process_data()).
    std::vector<std::pair<std::string, std::string>> familyAndName = {
        {"pq", "PQ 26/25"},
        {"etd", "ETD 29"},
        {"rm", "RM 8"},
        {"e", "E 25/13/7"},
        {"ec", "EC 35"},
        {"pq", "PQ27.3/18"},
    };

    for (auto& [family, name] : familyAndName) {
        DYNAMIC_SECTION("Standard shape reference: " << family << " / " << name) {
            json shapeJson = {{"type", "standard"}, {"family", family}, {"name", name}};
            json coreJson;
            coreJson["name"] = "t";
            coreJson["functionalDescription"] = {
                {"type", "pieceAndPlate"}, {"material", "3C95"}, {"shape", shapeJson},
                {"gapping", json::array()}, {"numberStacks", 1}};

            Core core(coreJson, false, false, false);
            REQUIRE_NOTHROW(core.process_data());

            REQUIRE(core.get_processed_description());
            auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
            CHECK(effectiveParameters.get_effective_area() > 0);
            CHECK(effectiveParameters.get_effective_length() > 0);
            CHECK(effectiveParameters.get_effective_volume() > 0);
        }
    }
}

TEST_CASE("Test_Core_Initial_Permeability", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);

    Core core(coreJson, true);

    auto initialPermeability = core.get_initial_permeability(25);
    double expectedInitialPermeability = 2270;

    REQUIRE_THAT(expectedInitialPermeability, Catch::Matchers::WithinAbs(initialPermeability, expectedInitialPermeability * maximumError));

    initialPermeability = core.get_initial_permeability(125);
    expectedInitialPermeability = 3975;
    REQUIRE_THAT(expectedInitialPermeability, Catch::Matchers::WithinAbs(initialPermeability, expectedInitialPermeability * maximumError));
}

TEST_CASE("Test_Core_Effective_Permeability", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);

    Core core(coreJson, true);

    auto effectivePermeability = core.get_effective_permeability(25);
    double expectedEffectivePermeability = 136;

    REQUIRE_THAT(expectedEffectivePermeability, Catch::Matchers::WithinAbs(effectivePermeability, expectedEffectivePermeability * maximumError));

    effectivePermeability = core.get_effective_permeability(125);
    expectedEffectivePermeability = 139;
    REQUIRE_THAT(expectedEffectivePermeability, Catch::Matchers::WithinAbs(effectivePermeability, expectedEffectivePermeability * maximumError));
}

TEST_CASE("Test_Core_Reluctance", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);

    Core core(coreJson, true);

    auto reluctance = core.get_reluctance(25);
    double expectedReluctance = 1.02e+06;

    REQUIRE_THAT(expectedReluctance, Catch::Matchers::WithinAbs(reluctance, expectedReluctance * maximumError));

    reluctance = core.get_reluctance(125);
    expectedReluctance = 997019;
    REQUIRE_THAT(expectedReluctance, Catch::Matchers::WithinAbs(reluctance, expectedReluctance * maximumError));
}

TEST_CASE("Test_Core_Resistivity", "[constructive-model][core][functional-description][smoke-test]") {
    auto coreFilePath = masPath + "samples/magnetic/core/core_E_55_21_N97_additive.json";
    std::ifstream json_file(coreFilePath);

    auto coreJson = json::parse(json_file);

    Core core(coreJson, true);

    auto resistivity = core.get_resistivity(25);
    double expectedResistivity = 8;

    REQUIRE_THAT(expectedResistivity, Catch::Matchers::WithinAbs(resistivity, expectedResistivity * maximumError));

    resistivity = core.get_resistivity(125);
    expectedResistivity = 8;
    REQUIRE_THAT(expectedResistivity, Catch::Matchers::WithinAbs(resistivity, expectedResistivity * maximumError));
}

TEST_CASE("Toroid_Effective_Parameters_Different_Standards", "[constructive-model][core][effective-parameters][standards]") {
    // Test that powder materials use IEC 63182 and ferrite materials use IEC 60205
    // and that the resulting effective parameters are different
    
    json coreJsonPowder;
    coreJsonPowder["functionalDescription"] = json();
    coreJsonPowder["name"] = "core_T_40_24_16_powder";
    coreJsonPowder["functionalDescription"]["type"] = "toroidal";
    coreJsonPowder["functionalDescription"]["material"] = "Kool Mµ 14";  // Powder material
    coreJsonPowder["functionalDescription"]["shape"] = "T 40/24/16";
    coreJsonPowder["functionalDescription"]["gapping"] = json::array();
    coreJsonPowder["functionalDescription"]["numberStacks"] = 1;
    
    json coreJsonFerrite;
    coreJsonFerrite["functionalDescription"] = json();
    coreJsonFerrite["name"] = "core_T_40_24_16_ferrite";
    coreJsonFerrite["functionalDescription"]["type"] = "toroidal";
    coreJsonFerrite["functionalDescription"]["material"] = "N97";  // Ferrite material
    coreJsonFerrite["functionalDescription"]["shape"] = "T 40/24/16";
    coreJsonFerrite["functionalDescription"]["gapping"] = json::array();
    coreJsonFerrite["functionalDescription"]["numberStacks"] = 1;
    
    // Create cores
    Core corePowder(coreJsonPowder, true);
    Core coreFerrite(coreJsonFerrite, true);
    
    // Verify material types
    REQUIRE(corePowder.resolve_material().get_material() == MAS::MaterialType::POWDER);
    REQUIRE(coreFerrite.resolve_material().get_material() == MAS::MaterialType::FERRITE);
    
    // Get effective parameters
    auto paramsPowder = corePowder.get_processed_description()->get_effective_parameters();
    auto paramsFerrite = coreFerrite.get_processed_description()->get_effective_parameters();
    
    // Verify that effective parameters are calculated
    REQUIRE(paramsPowder.get_effective_area() > 0);
    REQUIRE(paramsPowder.get_effective_length() > 0);
    REQUIRE(paramsPowder.get_effective_volume() > 0);
    REQUIRE(paramsPowder.get_minimum_area() > 0);
    
    REQUIRE(paramsFerrite.get_effective_area() > 0);
    REQUIRE(paramsFerrite.get_effective_length() > 0);
    REQUIRE(paramsFerrite.get_effective_volume() > 0);
    REQUIRE(paramsFerrite.get_minimum_area() > 0);
    
    // The key assertion: parameters should be DIFFERENT between standards
    // IEC 63182 (powder) should give different results than IEC 60205 (ferrite)
    bool areaDifferent = paramsPowder.get_effective_area() != paramsFerrite.get_effective_area();
    bool lengthDifferent = paramsPowder.get_effective_length() != paramsFerrite.get_effective_length();
    bool volumeDifferent = paramsPowder.get_effective_volume() != paramsFerrite.get_effective_volume();
    
    // At least one parameter should be different
    REQUIRE((areaDifferent || lengthDifferent || volumeDifferent));
}

TEST_CASE("Toroid_Coating_Winding_Window_Offset", "[constructive-model][core][coating][winding-window]") {
    // §5 geometry: a core coating lines the bore, so the usable winding bore
    // shrinks by the coating thickness, while the effective magnetic parameters
    // (computed from the bare ferrite) are unchanged. Uncoated cores are a no-op.
    auto makeToroid = [](const std::string& coating) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_T_40_24_16_coating";
        coreJson["functionalDescription"]["type"] = "toroidal";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "T 40/24/16";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        if (!coating.empty()) {
            coreJson["functionalDescription"]["coating"] = coating;
        }
        return Core(coreJson, true);
    };

    Core defaulted = makeToroid("");
    Core epoxy = makeToroid("epoxy");
    Core parylene = makeToroid("parylene");

    // A toroid with no coating field is jacketed in practice, so it falls back to the
    // default (epoxy) coating; a name-only coating resolves to its datasheet default.
    //
    // ABT #964: this is a FERRITE ring core, and epoxy on ferrite is not the film a powder
    // toroid carries. Fair-Rite's thermo-set plastic adds at most 0.5 mm across a diameter
    // (0.25 mm per surface), Ferroxcube DIMENSIONS its TN jacket at ~0.3 mm, and TDK's coated
    // limits on R 50.0x30.0x20.0 work out to exactly 0.400 mm per surface against the UNCOATED
    // limits, matching the "< 0.4 mm" it states. The old 0.10 mm here was the POWDER-core value,
    // which is still what a powder toroid resolves to.
    REQUIRE_THAT(defaulted.get_coating_thickness(),
                 Catch::Matchers::WithinAbs(Defaults().defaultFerriteEpoxyCoreCoatingThickness, 1e-12));
    REQUIRE_THAT(epoxy.get_coating_thickness(),
                 Catch::Matchers::WithinAbs(Defaults().defaultFerriteEpoxyCoreCoatingThickness, 1e-12));
    REQUIRE_THAT(parylene.get_coating_thickness(), Catch::Matchers::WithinAbs(12.7e-6, 1e-12));

    // True-uncoated baseline: an explicit zero-thickness coating gives the full bore.
    json uncoatedJson;
    uncoatedJson["functionalDescription"] = json();
    uncoatedJson["name"] = "core_T_40_24_16_uncoated";
    uncoatedJson["functionalDescription"]["type"] = "toroidal";
    uncoatedJson["functionalDescription"]["material"] = "N97";
    uncoatedJson["functionalDescription"]["shape"] = "T 40/24/16";
    uncoatedJson["functionalDescription"]["gapping"] = json::array();
    uncoatedJson["functionalDescription"]["numberStacks"] = 1;
    uncoatedJson["functionalDescription"]["coating"] = {{"type", "epoxy"}, {"thickness", 0.0}};
    Core uncoated(uncoatedJson, true);
    REQUIRE(uncoated.get_coating_thickness() == 0.0);

    double uncoatedRh = uncoated.get_processed_description()->get_winding_windows()[0].get_radial_height().value();
    double epoxyRh = epoxy.get_processed_description()->get_winding_windows()[0].get_radial_height().value();
    double paryleneRh = parylene.get_processed_description()->get_winding_windows()[0].get_radial_height().value();

    // The usable winding bore shrinks by exactly the coating thickness vs the bare bore.
    REQUIRE(epoxyRh < uncoatedRh);
    REQUIRE_THAT(uncoatedRh - epoxyRh,
                 Catch::Matchers::WithinAbs(Defaults().defaultFerriteEpoxyCoreCoatingThickness, 1e-9));
    REQUIRE_THAT(uncoatedRh - paryleneRh, Catch::Matchers::WithinAbs(12.7e-6, 1e-9));
    // Window area follows the reduced radius.
    double epoxyArea = epoxy.get_processed_description()->get_winding_windows()[0].get_area().value();
    REQUIRE_THAT(epoxyArea, Catch::Matchers::WithinRel(std::numbers::pi * epoxyRh * epoxyRh, 1e-9));

    // Effective magnetic parameters are bare-ferrite and must be UNCHANGED by coating.
    auto pb = uncoated.get_processed_description()->get_effective_parameters();
    auto pe = epoxy.get_processed_description()->get_effective_parameters();
    REQUIRE(pb.get_effective_area() == pe.get_effective_area());
    REQUIRE(pb.get_effective_length() == pe.get_effective_length());
    REQUIRE(pb.get_effective_volume() == pe.get_effective_volume());

    // An explicit coating thickness (object form) overrides the name default.
    json explicitJson;
    explicitJson["functionalDescription"] = json();
    explicitJson["name"] = "core_T_40_24_16_explicit_coating";
    explicitJson["functionalDescription"]["type"] = "toroidal";
    explicitJson["functionalDescription"]["material"] = "N97";
    explicitJson["functionalDescription"]["shape"] = "T 40/24/16";
    explicitJson["functionalDescription"]["gapping"] = json::array();
    explicitJson["functionalDescription"]["numberStacks"] = 1;
    explicitJson["functionalDescription"]["coating"] = {{"type", "epoxy"}, {"thickness", 0.0005}};
    Core explicitCore(explicitJson, true);
    REQUIRE_THAT(explicitCore.get_coating_thickness(), Catch::Matchers::WithinAbs(0.0005, 1e-12));
    double explicitRh = explicitCore.get_processed_description()->get_winding_windows()[0].get_radial_height().value();
    REQUIRE_THAT(uncoatedRh - explicitRh, Catch::Matchers::WithinAbs(0.0005, 1e-9));

    // A NON-toroidal core is wound on a bobbin, not directly on the ferrite, so it gets
    // NO default coating even when the coating field is absent.
    json eCoreJson;
    eCoreJson["functionalDescription"] = json();
    eCoreJson["name"] = "core_E_no_coating";
    eCoreJson["functionalDescription"]["type"] = "two-piece set";
    eCoreJson["functionalDescription"]["material"] = "N97";
    eCoreJson["functionalDescription"]["shape"] = "E 42/21/20";
    eCoreJson["functionalDescription"]["gapping"] = json::array();
    eCoreJson["functionalDescription"]["numberStacks"] = 1;
    Core eCore(eCoreJson, true);
    REQUIRE(eCore.get_coating_thickness() == 0.0);

    // A SMALL toroid (OD <= 5.08 mm) defaults to parylene rather than epoxy, mirroring
    // manufacturer practice (thin film on small bores). The large T 40/24/16 above gets
    // epoxy (0.1 mm); this small one gets parylene (12.7 um) with no coating field.
    json smallJson;
    smallJson["functionalDescription"] = json();
    smallJson["name"] = "core_T_small_parylene_default";
    smallJson["functionalDescription"]["type"] = "toroidal";
    smallJson["functionalDescription"]["material"] = "N97";
    smallJson["functionalDescription"]["shape"] = "T 2.5/1.5/1";
    smallJson["functionalDescription"]["gapping"] = json::array();
    smallJson["functionalDescription"]["numberStacks"] = 1;
    Core smallToroid(smallJson, true);
    REQUIRE_THAT(smallToroid.get_coating_thickness(), Catch::Matchers::WithinAbs(12.7e-6, 1e-12));
    REQUIRE_THAT(smallToroid.get_coating_relative_permittivity(), Catch::Matchers::WithinAbs(3.1, 1e-9));
}

TEST_CASE("Toroid_Coating_Relative_Permittivity", "[constructive-model][core][coating]") {
    // §7 prerequisite: the coating's relative permittivity is resolved from its
    // insulation material — datasheet default per type, or an explicit material.
    auto makeCore = [](const json& coating) {
        json coreJson;
        coreJson["functionalDescription"] = json();
        coreJson["name"] = "core_coating_permittivity";
        coreJson["functionalDescription"]["type"] = "toroidal";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = "T 40/24/16";
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        if (!coating.is_null()) {
            coreJson["functionalDescription"]["coating"] = coating;
        }
        return Core(coreJson, true);
    };

    // Name-only coatings resolve to the datasheet default material permittivity.
    REQUIRE_THAT(makeCore("epoxy").get_coating_relative_permittivity(), Catch::Matchers::WithinAbs(3.6, 1e-9));
    REQUIRE_THAT(makeCore("parylene").get_coating_relative_permittivity(), Catch::Matchers::WithinAbs(3.1, 1e-9));

    // An explicit material overrides the per-type default (Kapton HN = 3.4).
    Core explicitMat = makeCore(json({{"type", "epoxy"}, {"thickness", 0.0003}, {"material", "Kapton HN"}}));
    REQUIRE_THAT(explicitMat.get_coating_relative_permittivity(), Catch::Matchers::WithinAbs(3.4, 1e-9));

    // A toroid with no coating field falls back to the default epoxy permittivity.
    REQUIRE_THAT(makeCore(json(nullptr)).get_coating_relative_permittivity(), Catch::Matchers::WithinAbs(3.6, 1e-9));

    // A NON-toroidal core has no coating, so resolving its permittivity throws rather
    // than inventing one (it is wound on a bobbin, not directly on the ferrite).
    json eCoreJson;
    eCoreJson["functionalDescription"] = json();
    eCoreJson["name"] = "core_E_no_coating_permittivity";
    eCoreJson["functionalDescription"]["type"] = "two-piece set";
    eCoreJson["functionalDescription"]["material"] = "N97";
    eCoreJson["functionalDescription"]["shape"] = "E 42/21/20";
    eCoreJson["functionalDescription"]["gapping"] = json::array();
    eCoreJson["functionalDescription"]["numberStacks"] = 1;
    REQUIRE_THROWS(Core(eCoreJson, true).get_coating_relative_permittivity());
}

}  // namespace

namespace TestDrumCore {
    // ABT #331: drum (open-shape) core geometry. The winding window is the groove between the
    // flanges, the single column is the post with the bore subtracted, and nothing is doubled.
    TEST_CASE("Test_Drum_Core_Geometry", "[core][drum][open-core]") {
        settings.reset();
        clear_databases();
        auto core = OpenMagneticsTesting::get_quick_core("DRH-14X20-4C", json::array(), 1, "Dummy");
        REQUIRE(core.get_functional_description().get_type() == CoreType::OPEN_SHAPE);
        auto processed = core.get_processed_description().value();

        auto windingWindow = processed.get_winding_windows()[0];
        // groove: width (A-C)/2 = (14-9)/2 = 2.5 mm, height E = 12.5 mm
        CHECK_THAT(windingWindow.get_width().value(), Catch::Matchers::WithinRel(0.0025, 1e-6));
        CHECK_THAT(windingWindow.get_height().value(), Catch::Matchers::WithinRel(0.0125, 1e-6));

        auto columns = processed.get_columns();
        REQUIRE(columns.size() == 1);
        // post area pi/4 (C^2 - H^2) = pi/4 (9^2 - 3.2^2) mm^2
        double expectedArea = std::numbers::pi / 4 * (pow(0.009, 2) - pow(0.0032, 2));
        CHECK_THAT(columns[0].get_area(), Catch::Matchers::WithinRel(expectedArea, 0.01));
        // height NOT doubled: single piece
        CHECK_THAT(processed.get_height(), Catch::Matchers::WithinRel(0.020, 1e-6));
        settings.reset();
    }

    // ABT #366: shielded drum (drumRing) — a drum closed by a concentric shield ring,
    // CoreType::PIECE_AND_PLATE. The winding window stays the drum groove, the envelope
    // includes the ring, and process_gap synthesizes the two STRUCTURAL annular clearance
    // gaps ((K - A)/2 each, unrolled-annulus section) instead of the column-based machinery.
    TEST_CASE("Test_Drum_Ring_Core_Geometry", "[core][drum-ring]") {
        settings.reset();
        clear_databases();
        auto core = OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "Dummy");
        REQUIRE(core.get_functional_description().get_type() == CoreType::PIECE_AND_PLATE);
        auto processed = core.get_processed_description().value();

        // Winding window = drum groove: width (A - C)/2 = (2.3 - 1.1)/2 mm, height E = 0.58 mm.
        auto windingWindow = processed.get_winding_windows()[0];
        CHECK_THAT(windingWindow.get_width().value(), Catch::Matchers::WithinRel((0.0023 - 0.0011) / 2, 1e-6));
        CHECK_THAT(windingWindow.get_height().value(), Catch::Matchers::WithinRel(0.00058, 1e-6));

        // Envelope includes the ring: width/depth = ring OD J, height = max(B, L).
        CHECK_THAT(processed.get_width(), Catch::Matchers::WithinRel(0.003, 1e-6));
        CHECK_THAT(processed.get_height(), Catch::Matchers::WithinRel(0.00105, 1e-6));

        // Two structural annular gaps: residual, length (K - A)/2 = 50 um, area = mean
        // cylindrical surface over each flange thickness, mirrored about the equator.
        auto gapping = core.get_functional_description().get_gapping();
        REQUIRE(gapping.size() == 2);
        double meanRadius = (0.0023 + 0.0024) / 4;
        for (auto& gap : gapping) {
            CHECK(gap.get_type() == GapType::RESIDUAL);
            CHECK_THAT(gap.get_length(), Catch::Matchers::WithinRel((0.0024 - 0.0023) / 2, 1e-4));
            CHECK_THAT(gap.get_area().value(),
                       Catch::Matchers::WithinRel(2 * std::numbers::pi * meanRadius * 0.00021, 0.01));
        }
        CHECK_THAT(gapping[0].get_coordinates().value()[1],
                   Catch::Matchers::WithinAbs((0.001 - 0.00021) / 2, 1e-6));
        CHECK_THAT(gapping[1].get_coordinates().value()[1],
                   Catch::Matchers::WithinAbs(-(0.001 - 0.00021) / 2, 1e-6));

        // The geometrical description is drum solid (CLOSED) + ring closer (PLATE).
        auto geometricalDescription = core.create_geometrical_description().value();
        REQUIRE(geometricalDescription.size() == 2);
        CHECK(geometricalDescription[0].get_type() == CoreGeometricalDescriptionElementType::CLOSED);
        CHECK(geometricalDescription[1].get_type() == CoreGeometricalDescriptionElementType::PLATE);

        // User gapping on a drumRing is rejected: nothing can be ground on the assembly.
        auto userGapping = json::array();
        userGapping.push_back(json{{"type", "subtractive"}, {"length", 0.0001}});
        CHECK_THROWS(OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", userGapping, 1, "Dummy"));
        settings.reset();
    }

    // ABT #357: molded composite body (WE-MAPI class) — a single pressed CLOSED solid whose
    // distributed gap lives in the material. Magnetically a pot core with a rectangular outer
    // boundary: post + two plates + return shell. MAPI-4020-like custom dimensions (vendors
    // publish no internals; the cavity here is a plausible reconstruction for geometry tests).
    // Letters per the pot-core convention: D cavity height, E cavity OD, F post diameter.
    TEST_CASE("Test_Molded_Core_Geometry", "[core][molded]") {
        settings.reset();
        clear_databases();
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()}, {"name", "MAPI-like 4020"},
            {"dimensions", {
                {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
                {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}
        };
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", "Kool Mµ 26"}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        Core core(coreJson);
        core.process_data();
        core.process_gap();
        REQUIRE(core.get_functional_description().get_type() == CoreType::CLOSED_SHAPE);
        auto processed = core.get_processed_description().value();

        // Winding window = the coil cavity annulus (pot-core letter convention):
        // width (E - F)/2, height D.
        auto windingWindow = processed.get_winding_windows()[0];
        CHECK_THAT(windingWindow.get_width().value(), Catch::Matchers::WithinRel((0.0030 - 0.0012) / 2, 1e-6));
        CHECK_THAT(windingWindow.get_height().value(), Catch::Matchers::WithinRel(0.0014, 1e-6));

        // Central post + return-shell columns.
        auto columns = processed.get_columns();
        REQUIRE(columns.size() == 2);
        CHECK_THAT(columns[0].get_area(),
                   Catch::Matchers::WithinRel(std::numbers::pi / 4 * pow(0.0012, 2), 0.01));
        double expectedShellArea = 0.0041 * 0.0041 - std::numbers::pi / 4 * pow(0.0030, 2);
        CHECK_THAT(columns[1].get_area(), Catch::Matchers::WithinRel(expectedShellArea, 0.01));

        // Body envelope, nothing doubled.
        CHECK_THAT(processed.get_width(), Catch::Matchers::WithinRel(0.0041, 1e-6));
        CHECK_THAT(processed.get_height(), Catch::Matchers::WithinRel(0.0021, 1e-6));

        // No gaps, ever: neither synthesized residual gaps nor user gapping.
        CHECK(core.get_functional_description().get_gapping().empty());
        auto effectiveParameters = processed.get_effective_parameters();
        double bodyVolume = 0.0041 * 0.0041 * 0.0021;
        CHECK(effectiveParameters.get_effective_volume() > 0);
        CHECK(effectiveParameters.get_effective_volume() < bodyVolume);
        // The path length must exceed one cavity-height-plus-radial loop and stay below a few
        // body perimeters — a gross sectioning error would leave this band.
        CHECK(effectiveParameters.get_effective_length() > 0.0021);
        CHECK(effectiveParameters.get_effective_length() < 4 * (0.0041 + 0.0021));

        // Single CLOSED solid in the geometrical description.
        auto geometricalDescription = core.create_geometrical_description().value();
        REQUIRE(geometricalDescription.size() == 1);
        CHECK(geometricalDescription[0].get_type() == CoreGeometricalDescriptionElementType::CLOSED);

        // Discrete gapping is physically meaningless on a molded body. The Core(json)
        // constructor already processes, so the throw fires during construction.
        json gappedCoreJson = coreJson;
        gappedCoreJson["functionalDescription"]["gapping"].push_back(
            json{{"type", "subtractive"}, {"length", 0.0001}});
        CHECK_THROWS(Core(gappedCoreJson).process_gap());
        settings.reset();
    }

    // ABT #362: semi-shielded drum — a wound drum overcoated with magnetic epoxy acting as a
    // low-mu return shell (WE-LQS class). Letters: drum A..H + shell envelope J/K/L. The shell
    // is cast in contact (no gaps); its material rides the magneticEpoxy coating.
    TEST_CASE("Test_Drum_Semishielded_Core_Geometry", "[core][drum-semishielded]") {
        settings.reset();
        clear_databases();
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
            {"aliases", json::array()}, {"name", "LQS-like 4018"},
            {"dimensions", {
                {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
                {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}},
                {"J", {{"nominal", 0.0040}}}, {"K", {{"nominal", 0.0040}}}, {"L", {{"nominal", 0.0018}}}}}
        };
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "pieceAndPlate"}, {"material", "3C90"}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1},
            {"coating", {{"type", "magneticEpoxy"}, {"thickness", 0.0001}, {"material", "Kool Mµ 26"}}}};
        Core core(coreJson);
        core.process_data();
        core.process_gap();
        REQUIRE(core.get_functional_description().get_type() == CoreType::PIECE_AND_PLATE);
        auto processed = core.get_processed_description().value();

        // Winding window stays the drum groove: width (A - C)/2, height E.
        auto windingWindow = processed.get_winding_windows()[0];
        CHECK_THAT(windingWindow.get_width().value(), Catch::Matchers::WithinRel((0.0038 - 0.0015) / 2, 1e-6));
        CHECK_THAT(windingWindow.get_height().value(), Catch::Matchers::WithinRel(0.0010, 1e-6));

        // Envelope = the finished body J x K x L.
        CHECK_THAT(processed.get_width(), Catch::Matchers::WithinRel(0.0040, 1e-6));
        CHECK_THAT(processed.get_height(), Catch::Matchers::WithinRel(0.0018, 1e-6));

        // No gaps ever: the glue is cast in contact.
        CHECK(core.get_functional_description().get_gapping().empty());
        auto effectiveParameters = processed.get_effective_parameters();
        CHECK(effectiveParameters.get_effective_length() > 0.0018);
        CHECK(effectiveParameters.get_effective_length() < 4 * (0.0040 + 0.0018));

        // Drum solid (CLOSED) + glue shell closer (PLATE).
        auto geometricalDescription = core.create_geometrical_description().value();
        REQUIRE(geometricalDescription.size() == 2);
        CHECK(geometricalDescription[0].get_type() == CoreGeometricalDescriptionElementType::CLOSED);
        CHECK(geometricalDescription[1].get_type() == CoreGeometricalDescriptionElementType::PLATE);

        // User gapping rejected.
        json gappedCoreJson = coreJson;
        gappedCoreJson["functionalDescription"]["gapping"].push_back(
            json{{"type", "subtractive"}, {"length", 0.0001}});
        CHECK_THROWS(Core(gappedCoreJson).process_gap());

        // A shell envelope smaller than the drum is impossible (letters describe the FINISHED body).
        json badShapeCoreJson = coreJson;
        badShapeCoreJson["functionalDescription"]["shape"]["dimensions"]["J"] = {{"nominal", 0.0030}};
        CHECK_THROWS(Core(badShapeCoreJson).process_data());
        settings.reset();
    }

    // ABT #366/#362/#357: the 2D painter used to reconstruct a mirrored half-set for EVERY
    // non-toroidal family, drawing drums/drumRings/semi-shielded/molded as nonsense; they now
    // route to paint_drum_family_core. Smoke: each family paints a valid, non-empty SVG.
    TEST_CASE("Test_Drum_Family_Core_Painter_Smoke", "[core][drum][drum-ring][drum-semishielded][molded][painter]") {
        settings.reset();
        clear_databases();
        auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}
                                  .parent_path().append("..").append("output");
        std::filesystem::create_directories(outputFilePath);

        json semishieldedShape = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
            {"aliases", json::array()}, {"name", "LQS-like 4018"},
            {"dimensions", {
                {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
                {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}},
                {"J", {{"nominal", 0.0040}}}, {"K", {{"nominal", 0.0040}}}, {"L", {{"nominal", 0.0018}}}}}
        };
        json moldedShape = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()}, {"name", "MAPI-like 4020"},
            {"dimensions", {
                {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
                {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}
        };

        std::vector<std::pair<std::string, Core>> cores;
        cores.emplace_back("drum", OpenMagneticsTesting::get_quick_core("DRH-14X20-4C", json::array(), 1, "Dummy"));
        cores.emplace_back("drum_ring", OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "Dummy"));
        for (auto& [label, shapeJson] : std::vector<std::pair<std::string, json>>{
                 {"drum_semishielded", semishieldedShape}, {"molded", moldedShape}}) {
            json coreJson;
            coreJson["functionalDescription"] = {
                {"type", label == "molded" ? "closedShape" : "pieceAndPlate"}, {"material", "Dummy"},
                {"shape", shapeJson}, {"gapping", json::array()}, {"numberStacks", 1}};
            Core core(coreJson);
            core.process_data();
            cores.emplace_back(label, core);
        }

        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array();
        coilJson["functionalDescription"].push_back(json{
            {"name", "winding 0"}, {"numberTurns", 1}, {"numberParallels", 1},
            {"isolationSide", "primary"}, {"wire", "Dummy"}});

        for (auto& [label, core] : cores) {
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_" + label + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_core(magnetic);
            painter.export_svg();
            OpenMagneticsTesting::check_svg(outFile);
        }
        settings.reset();
    }
}

// ABT #366/#362/#357: SATURATION CURRENT for the new families. Isat is a headline datasheet
// spec for every one of them (drum, shielded drum, semi-shielded, molded), and it runs through
// a different chain than inductance — Bsat(T) x N x Ae / L — where L itself comes from the
// family's own model (open-core demagnetising for a bare drum, mixed-material sectioning for a
// semi-shielded, closed-circuit for drumRing/molded). None of that was exercised, so this pins
// that the chain runs and orders physically: for one turn count, MORE inductance means LESS
// saturation current, so the shielded assembly must saturate EARLIER than the bare drum it is
// built from.
namespace TestNewFamilySaturation {
    TEST_CASE("Test_Drum_Family_Saturation_Current", "[core][drum][drum-ring][molded][saturation]") {
        settings.reset();
        clear_databases();
        int64_t numberTurns = 20;

        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array({{
            {"name", "winding 0"}, {"numberTurns", numberTurns}, {"numberParallels", 1},
            {"isolationSide", "primary"}, {"wire", "Dummy"}}});

        auto magneticFromCore = [&](const Core& core) {
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(OpenMagnetics::Coil(coilJson));
            return magnetic;
        };

        // Bare drum with the same drum dimensions as the shielded pair below.
        json bareShape = {
            {"magneticCircuit", "open"}, {"type", "custom"}, {"family", "drum"},
            {"aliases", json::array()}, {"name", "DR 2.3 bare"},
            {"dimensions", {
                {"A", {{"nominal", 0.0023}}}, {"B", {{"nominal", 0.001}}}, {"C", {{"nominal", 0.0011}}},
                {"D", {{"nominal", 0.00021}}}, {"E", {{"nominal", 0.00058}}}, {"F", {{"nominal", 0.00021}}}}}
        };
        json bareCoreJson;
        bareCoreJson["functionalDescription"] = {
            {"type", "openShape"}, {"material", "3C90"}, {"shape", bareShape},
            {"gapping", json::array()}, {"numberStacks", 1}};
        Core bareDrum(bareCoreJson);
        bareDrum.process_data();
        auto bareMagnetic = magneticFromCore(bareDrum);
        double bareSaturationCurrent = bareMagnetic.calculate_saturation_current(25);

        // Shielded drum: same drum, closed by its ring.
        auto shieldedCore = OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90");
        auto shieldedMagnetic = magneticFromCore(shieldedCore);
        double shieldedSaturationCurrent = shieldedMagnetic.calculate_saturation_current(25);

        // Molded composite body.
        json moldedShape = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()}, {"name", "MAPI-like 4020"},
            {"dimensions", {
                {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
                {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}
        };
        json moldedCoreJson;
        moldedCoreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", "Kool Mµ 26"}, {"shape", moldedShape},
            {"gapping", json::array()}, {"numberStacks", 1}};
        Core moldedCore(moldedCoreJson);
        moldedCore.process_data();
        moldedCore.process_gap();
        auto moldedMagnetic = magneticFromCore(moldedCore);
        double moldedSaturationCurrent = moldedMagnetic.calculate_saturation_current(25);

        // Semi-shielded: same drum, closed by a magnetic-epoxy shell (its permeability comes
        // from the magneticEpoxy coating, so this also proves the mixed-material inductance path
        // feeds the saturation chain).
        json semishieldedShape = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
            {"aliases", json::array()}, {"name", "LQS-like 4018"},
            {"dimensions", {
                {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
                {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}},
                {"J", {{"nominal", 0.0040}}}, {"K", {{"nominal", 0.0040}}}, {"L", {{"nominal", 0.0018}}}}}
        };
        json semishieldedCoreJson;
        semishieldedCoreJson["functionalDescription"] = {
            {"type", "pieceAndPlate"}, {"material", "3C90"}, {"shape", semishieldedShape},
            {"gapping", json::array()}, {"numberStacks", 1},
            {"coating", {{"type", "magneticEpoxy"}, {"thickness", 0.0001}, {"material", "Kool Mµ 26"}}}};
        Core semishieldedCore(semishieldedCoreJson);
        semishieldedCore.process_data();
        semishieldedCore.process_gap();
        auto semishieldedMagnetic = magneticFromCore(semishieldedCore);
        double semishieldedSaturationCurrent = semishieldedMagnetic.calculate_saturation_current(25);

        UNSCOPED_INFO("Isat @20T: bare drum " << bareSaturationCurrent << " A, shielded "
                      << shieldedSaturationCurrent << " A, semi-shielded "
                      << semishieldedSaturationCurrent << " A, molded " << moldedSaturationCurrent << " A");
        for (double saturationCurrent : {bareSaturationCurrent, shieldedSaturationCurrent,
                                         semishieldedSaturationCurrent, moldedSaturationCurrent}) {
            CHECK(std::isfinite(saturationCurrent));
            CHECK(saturationCurrent > 0);
        }
        // The ferrite ring multiplies inductance (~2.9x on this pair), so with the same turns the
        // shielded assembly reaches Bsat at a proportionally lower current than the bare drum.
        CHECK(shieldedSaturationCurrent < bareSaturationCurrent);
        settings.reset();
    }
}

// ABT #379: a caller that DECLARES a multicolumn core (one winding window per column) must keep
// that topology through processing. It used to be silently reduced to the single window the E
// family rebuilds from its functionalDescription, and the failure then surfaced far away and
// blamed the wrong object: the coil's section-derived placement resolved "window 2",
// ReluctanceNetwork found one window, and the error named the WINDING
// ("Winding Secondary references winding window 2 but the core has 1 winding windows") for
// something the CORE processing had dropped. Fixture mirrors OMFEM's multicolumn E42.
namespace TestMulticolumnWindows {
    TEST_CASE("Test_Core_Declared_Multicolumn_Windows_Survive_Processing", "[core][multicolumn]") {
        settings.reset();
        clear_databases();

        auto declaredCore = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::array(), 1, "3C97");
        settings.set_core_per_column_winding_windows(true);
        auto multicolumnCore = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::array(), 1, "3C97");
        auto declaredWindows = multicolumnCore.get_processed_description().value().get_winding_windows();
        settings.reset();
        REQUIRE(declaredWindows.size() > 1);  // the per-column machinery produced the topology

        // Hand that multi-window description to a core built the ordinary (single-window) way,
        // exactly as a caller supplying a multicolumn MAS file does, then re-process.
        auto processedDescription = declaredCore.get_processed_description().value();
        processedDescription.set_winding_windows(declaredWindows);
        declaredCore.set_processed_description(processedDescription);
        REQUIRE(declaredCore.get_processed_description().value().get_winding_windows().size() == declaredWindows.size());

        declaredCore.process_data();

        auto survivingWindows = declaredCore.get_processed_description().value().get_winding_windows();
        UNSCOPED_INFO("declared " << declaredWindows.size() << " winding windows, kept "
                      << survivingWindows.size() << " after processing");
        CHECK(survivingWindows.size() == declaredWindows.size());
        // Each window still names the column it wraps, which is what winding placement resolves.
        for (auto& windingWindow : survivingWindows) {
            CHECK(windingWindow.get_column().has_value());
        }
        settings.reset();
    }

    // The end-to-end shape of ABT #379, with the MAS file that reported it: a hand-authored
    // 3-column E 42 transformer whose windings carry no explicit windingWindow (placement lives
    // on the sections). Before the fix, re-processing collapsed 3 windows to 1, the coil's
    // section-derived window 2 then fell outside the core, and the thrown message blamed the
    // winding. Now the declared topology survives and the magnetic simulates.
    TEST_CASE("Test_Core_Multicolumn_Mas_File_Simulates", "[core][multicolumn]") {
        settings.reset();
        auto path = std::filesystem::path{std::source_location::current().file_name()}
                        .parent_path().append("testData").append("multicolumn_e42_transformer.json");
        std::ifstream masFile(path);
        REQUIRE(masFile.good());
        json masJson = json::parse(masFile);

        OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
        auto processedWindows = magnetic.get_core().get_processed_description().value().get_winding_windows();
        UNSCOPED_INFO("core kept " << processedWindows.size() << " winding windows after construction");
        CHECK(processedWindows.size() == 3);

        // The path that used to throw with the misleading message.
        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        OpenMagnetics::MagneticSimulator simulator;
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = simulator.simulate(inputs, magnetic));
        REQUIRE(simulated.get_outputs().size() > 0);

        // ABT #925: this fixture shipped with designRequirements.turnsRatios = 1e9 for a 24:12
        // coil, so its generated secondary excitation carried 1e9 A. Nothing noticed until
        // simulate() grew a thermal step (ABT #906): 1.26e16 W of ohmic loss in the secondary
        // drove the thermal network to 5.7e6 K and it refused to converge, and the throw read as
        // a thermal-solver bug rather than as impossible input. Pin the ampere-turn balance so a
        // regenerated fixture cannot go unphysical again without saying so here.
        auto primaryExcitation = inputs.get_winding_excitation(0, 0);
        auto secondaryExcitation = inputs.get_winding_excitation(0, 1);
        double primaryAmpereTurns = primaryExcitation.get_current()->get_processed()->get_rms().value()
                                    * magnetic.get_coil().get_functional_description()[0].get_number_turns();
        double secondaryAmpereTurns = secondaryExcitation.get_current()->get_processed()->get_rms().value()
                                      * magnetic.get_coil().get_functional_description()[1].get_number_turns();
        UNSCOPED_INFO("primary " << primaryAmpereTurns << " At, secondary " << secondaryAmpereTurns << " At");
        CHECK_THAT(secondaryAmpereTurns, Catch::Matchers::WithinRel(primaryAmpereTurns, 0.01));
        auto& output = simulated.get_outputs()[0];
        REQUIRE(output.get_inductance());
        double inductance = OpenMagnetics::resolve_dimensional_values(
            output.get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance());
        UNSCOPED_INFO("magnetizing inductance " << inductance * 1e6 << " uH");
        CHECK(std::isfinite(inductance));
        CHECK(inductance > 0);

        // Preserving the windows is only half the point: the network must actually USE them.
        // The primary sits on a LATERAL column, so its driving-point reluctance is
        // R_lateral + (R_central || R_other_lateral) -- strictly worse than the central-column
        // path the lumped N^2/R model assumes. Pin that the placed answer is the lower one, so a
        // future regression that silently reverts to the ideal single-window circuit is caught
        // even though nothing throws.
        // Preserving the windows is only half the point: the placement must actually REACH the
        // physics. The two windings sit on DIFFERENT legs of the E core, so they must resolve to
        // different columns -- with the collapsed single window this resolution is exactly what
        // threw. Their coupling is then the real flux divider (secondary flux splits between the
        // centre leg and the far leg), not the ideal rank-1 coupling a one-window core implies.
        REQUIRE(OpenMagnetics::ReluctanceNetwork::has_non_main_placement(magnetic));
        auto columnIndexPerWinding = OpenMagnetics::ReluctanceNetwork::resolve_winding_column_indexes(magnetic);
        REQUIRE(columnIndexPerWinding.size() == 2);
        UNSCOPED_INFO("Primary on column " << columnIndexPerWinding[0] << ", Secondary on column "
                      << columnIndexPerWinding[1]);
        CHECK(columnIndexPerWinding[0] != columnIndexPerWinding[1]);

        OpenMagnetics::Inductance inductanceModel;
        auto inductanceMatrix = inductanceModel.calculate_inductance_matrix(magnetic, 100000).get_magnitude();
        double selfPrimary = inductanceMatrix["Primary"]["Primary"].get_nominal().value();
        double selfSecondary = inductanceMatrix["Secondary"]["Secondary"].get_nominal().value();
        double mutual = inductanceMatrix["Primary"]["Secondary"].get_nominal().value();
        double couplingCoefficient = std::abs(mutual) / std::sqrt(selfPrimary * selfSecondary);
        UNSCOPED_INFO("coupling coefficient " << couplingCoefficient);
        CHECK(couplingCoefficient > 0.1);
        CHECK(couplingCoefficient < 0.999);
        settings.reset();
    }
}

// A CATALOGUE record whose gapping cannot fit its columns must be refused at load, naming itself.
// process_gap() reports that by returning false, which is a NORMAL answer during core advising —
// the CoreAdviser sweeps candidate gap lengths and necessarily generates some that do not fit,
// then skips them — but never legitimate for a shipped part. Left unchecked such a record carried
// gaps with no area and killed the first consumer that swept the whole catalogue, with a bare
// "[GAP_INVALID_DIMENSIONS] Gap Area is not set" naming neither the core nor the reason. Finding
// the culprits took a full-catalogue scan: seven Magnetics parts whose gap lengths were exact mil
// values stored 1000x too large, e.g. a 127 mm gap on a 19.3 mm core (5 mil, i.e. 0.127 mm). Data
// fixed in MAS; this keeps the next such record from getting in quietly.
TEST_CASE("Test_Catalogue_Core_With_Impossible_Gapping_Is_Refused_By_Name", "[core][gapping]") {
    settings.reset();
    clear_databases();

    // One record, shaped exactly like the catalogue bug: 127 mm of gap in a 19.3 mm core.
    json impossibleRecord;
    impossibleRecord["name"] = "E 19.3/4.8 - 3C97 - Gapped 127.000 mm";
    impossibleRecord["functionalDescription"] = json();
    impossibleRecord["functionalDescription"]["type"] = "twoPieceSet";
    impossibleRecord["functionalDescription"]["material"] = "3C97";
    impossibleRecord["functionalDescription"]["shape"] = "E 19.3/4.8";
    impossibleRecord["functionalDescription"]["numberStacks"] = 1;
    impossibleRecord["functionalDescription"]["gapping"] = json::array({
        {{"type", "subtractive"}, {"length", 0.127}},
        {{"type", "residual"}, {"length", 0.000005}},
        {{"type", "residual"}, {"length", 0.000005}},
    });

    std::string message;
    try {
        load_cores(impossibleRecord.dump());
        message = "no exception";
    }
    catch (const std::exception& exception) {
        message = exception.what();
    }
    UNSCOPED_INFO(message);
    CHECK(message.find("E 19.3/4.8 - 3C97 - Gapped 127.000 mm") != std::string::npos);  // names itself
    CHECK(message.find("does not fit its columns") != std::string::npos);

    // The same record with the gap it was meant to have loads, and its gaps carry an area.
    clear_databases();
    impossibleRecord["name"] = "E 19.3/4.8 - 3C97 - Gapped 0.127 mm";
    impossibleRecord["functionalDescription"]["gapping"][0]["length"] = 0.000127;
    REQUIRE_NOTHROW(load_cores(impossibleRecord.dump()));
    REQUIRE(OpenMagnetics::coreDatabase.size() == 1);
    for (auto& gap : OpenMagnetics::coreDatabase[0].get_functional_description().get_gapping()) {
        CHECK(gap.get_area().has_value());
    }
    clear_databases();
    settings.reset();
}

// ABT #680: calculate_core_gapping() and every other caller that hands process_gap() one
// specific core it expects to already be valid must not get back a schema-invalid gap (every
// derived field null) when the gap is too long for its column -- it must throw right where the
// mismatch is known. process_gap() itself must keep returning bare false: CoreAdviser sweeps
// candidate gap lengths and depends on that bool to silently skip the ones that don't fit
// (see the comment on the Core(json) constructor). process_gap_or_throw() is the loud variant.
TEST_CASE("ABT680_Gap_Longer_Than_Column_Opening_Throws_Instead_Of_Returning_Null_Fields", "[core][gapping]") {
    // Same custom E+I geometry and threshold as the bug report: the column opening is exactly
    // 2*D, confirmed there by bisection at several window heights.
    auto coreJson = [](double gapLength) {
        json j;
        j["name"] = "t";
        j["functionalDescription"]["type"] = "twoPieceSet";
        j["functionalDescription"]["material"] = "3C95";
        j["functionalDescription"]["shape"]["type"] = "custom";
        j["functionalDescription"]["shape"]["family"] = "ei";
        j["functionalDescription"]["shape"]["name"] = "t";
        j["functionalDescription"]["shape"]["dimensions"]["A"] = 0.0182;
        j["functionalDescription"]["shape"]["dimensions"]["B"] = 0.00745;
        j["functionalDescription"]["shape"]["dimensions"]["C"] = 0.0182;
        j["functionalDescription"]["shape"]["dimensions"]["D"] = 0.00047;
        j["functionalDescription"]["shape"]["dimensions"]["E"] = 0.01018;
        j["functionalDescription"]["shape"]["dimensions"]["F"] = 0.0014;
        j["functionalDescription"]["shape"]["dimensions"]["B2"] = 0.0014;
        j["functionalDescription"]["gapping"] = json::array({{{"type", "subtractive"}, {"length", gapLength}}});
        j["functionalDescription"]["numberStacks"] = 1;
        return j;
    };

    SECTION("just inside the threshold: process_gap() succeeds and every field is populated") {
        Core core(coreJson(0.00094), false, false, false);
        core.process_data();
        REQUIRE(core.process_gap());
        auto gapping = core.get_functional_description().get_gapping();
        REQUIRE(gapping.size() == 3u);
        CHECK(gapping[0].get_coordinates().has_value());
        CHECK(gapping[0].get_area().has_value());
    }

    SECTION("just past the threshold: process_gap() still just returns false") {
        Core core(coreJson(0.00095), false, false, false);
        core.process_data();
        CHECK_FALSE(core.process_gap());
    }

    SECTION("just past the threshold: process_gap_or_throw() throws, naming the gap and column") {
        Core core(coreJson(0.00095), false, false, false);
        core.process_data();
        std::string message = "no exception";
        try {
            core.process_gap_or_throw();
        }
        catch (const std::exception& exception) {
            message = exception.what();
        }
        UNSCOPED_INFO(message);
        CHECK(message != "no exception");
        CHECK(message.find("0.00095") != std::string::npos);
        CHECK(message.find("does not fit") != std::string::npos);
        // And it must not have left a schema-invalid gap (every derived field null) behind for
        // an inattentive caller to trip over downstream.
        auto gapping = core.get_functional_description().get_gapping();
        REQUIRE(gapping.size() == 1u);
        CHECK_FALSE(gapping[0].get_area().has_value());
    }
}

// Every core in the shipped catalogue must be constructible. This is what would have caught the
// seven bad records at the source instead of leaving them to break whichever consumer swept the
// full catalogue first (the core cross-referencer, which sets use_only_cores_in_stock(false)).
TEST_CASE("Test_All_Catalogue_Cores_Have_Feasible_Gapping", "[core][gapping][catalog]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    load_cores();
    auto& cores = OpenMagnetics::coreDatabase;
    REQUIRE(cores.size() > 1000);
    std::vector<std::string> coresWithoutGapArea;
    for (auto& core : cores) {
        for (auto& gap : core.get_functional_description().get_gapping()) {
            if (!gap.get_area()) {
                coresWithoutGapArea.push_back(core.get_name().value_or("<unnamed>"));
                break;
            }
        }
    }
    if (!coresWithoutGapArea.empty()) {
        std::string joined;
        for (size_t i = 0; i < coresWithoutGapArea.size() && i < 10; ++i) {
            joined += "\n    " + coresWithoutGapArea[i];
        }
        UNSCOPED_INFO(coresWithoutGapArea.size() << " catalogue cores have a gap with no area:" << joined);
    }
    CHECK(coresWithoutGapArea.empty());
    settings.reset();
}

// ABT #267/#407: a json -> Core conversion (json::parse(s).get<std::vector<Core>>(), or any
// nlohmann conversion) used to resolve through the base class to MAS's GENERATED
// from_json(json, MagneticCore&), skipping the legacy-form migration that the Core(json)
// constructor applies. The same document therefore produced different objects depending on which
// entry point you happened to use.
//
// It bit on a shape whose family was written "planar e", the pre-1.0 spelling of "planarE".
// Unmigrated, that string matches no enum value, and the generated converter leaves the enum
// DEFAULT-CONSTRUCTED — which is CoreShapeFamily::BLOCK, because block sorts first. So a planar E
// core silently became a "block" core, and the failure surfaced far away as "Unknown shape family:
// block" out of a factory that had never been asked for a block. (block is a declared-but-unused
// future family; nothing in the catalogue has one.)
TEST_CASE("Test_Core_Json_Conversion_Migrates_Legacy_Shape_Family", "[core][migration]") {
    settings.reset();
    clear_databases();

    json coreJson;
    coreJson["name"] = "legacy spelling core";
    coreJson["functionalDescription"] = json();
    coreJson["functionalDescription"]["type"] = "twoPieceSet";
    coreJson["functionalDescription"]["material"] = "3C97";
    coreJson["functionalDescription"]["numberStacks"] = 1;
    coreJson["functionalDescription"]["gapping"] = json::array();
    // The pre-1.0 spelling, exactly as it appears in older inventories.
    coreJson["functionalDescription"]["shape"] = json{
        {"name", "E 18/4/10"},
        {"family", "planar e"},
        {"type", "standard"},
        {"magneticCircuit", "open"},
        {"dimensions", json{
            {"A", {{"minimum", 0.01765}, {"maximum", 0.01835}}},
            {"B", {{"minimum", 0.0039},  {"maximum", 0.0041}}},
            {"C", {{"minimum", 0.0098},  {"maximum", 0.0102}}},
            {"D", {{"minimum", 0.0019},  {"maximum", 0.0021}}},
            {"E", {{"minimum", 0.0137},  {"maximum", 0.0143}}},
            {"F", {{"minimum", 0.0039},  {"maximum", 0.0041}}},
        }},
    };

    // Both entry points must agree, and both must say planar E rather than block.
    OpenMagnetics::Core constructedCore(coreJson);
    CHECK(constructedCore.get_shape_family() == CoreShapeFamily::PLANAR_E);

    auto convertedCore = coreJson.get<OpenMagnetics::Core>();
    CHECK(convertedCore.get_shape_family() == CoreShapeFamily::PLANAR_E);

    // And a family string that is not in the enum at all must be refused by name, rather than
    // quietly becoming value 0 the way "planar e" used to.
    json unknownFamilyJson = coreJson;
    unknownFamilyJson["functionalDescription"]["shape"]["family"] = "not a real family";
    std::string message;
    try {
        auto rejected = unknownFamilyJson.get<OpenMagnetics::Core>();
        message = "no exception";
    }
    catch (const std::exception& exception) {
        message = exception.what();
    }
    UNSCOPED_INFO(message);
    CHECK(message.find("not a real family") != std::string::npos);
    CHECK(message.find("E 18/4/10") != std::string::npos);
    settings.reset();
}

TEST_CASE("ABT644_Short_Gapping_List_Pads_With_Residual_Not_By_Repeating_Last",
          "[constructive-model][core][gapping]") {
    // A gapping list SHORTER than the column count used to be padded by repeating its LAST entry,
    // type included. One {subtractive, L} on a three-column core therefore became a subtractive gap
    // of length L in EVERY column -- a core ground on all three legs. The lateral gaps then sat in
    // parallel with the central one, so the reluctance was far too high and the core came back
    // 1.7-1.8x less inductive than the caller asked for, silently.
    //
    // The contract is that the gaps given map onto the columns in order (central column first) and
    // every remaining column gets a RESIDUAL gap. That is also exactly how all 2139 gapped cores in
    // the catalogue are written: [subtractive, residual, residual].
    auto constants = Constants();

    for (auto shapeName : {std::string("PQ 26/25"), std::string("E 42/21/20"),
                           std::string("ETD 29/16/10"), std::string("RM 10")}) {
        json coreJson;
        coreJson["name"] = "gapping test " + shapeName;
        coreJson["functionalDescription"] = json();
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "3C95";
        coreJson["functionalDescription"]["shape"] = shapeName;
        coreJson["functionalDescription"]["numberStacks"] = 1;
        coreJson["functionalDescription"]["gapping"] = json::array({
            json{{"type", "subtractive"}, {"length", 0.0005}}});

        OpenMagnetics::Core core(coreJson);
        auto gapping = core.get_functional_description().get_gapping();
        auto columns = core.get_processed_description().value().get_columns();

        UNSCOPED_INFO("shape: " << shapeName);
        // one gap per column, however many columns the shape has
        REQUIRE(gapping.size() == columns.size());
        REQUIRE(columns.size() > 1);

        // the stated gap lands on the central column...
        CHECK(columns[0].get_type() == ColumnType::CENTRAL);
        CHECK(gapping[0].get_type() == GapType::SUBTRACTIVE);
        CHECK_THAT(gapping[0].get_length(),
                   Catch::Matchers::WithinRel(0.0005, 1e-9));

        // ...and every other column gets a residual gap, NOT a copy of the subtractive one
        for (size_t i = 1; i < gapping.size(); ++i) {
            UNSCOPED_INFO("column index: " << i);
            CHECK(gapping[i].get_type() == GapType::RESIDUAL);
            CHECK_THAT(gapping[i].get_length(),
                       Catch::Matchers::WithinRel(constants.residualGap, 1e-9));
        }
    }
}

TEST_CASE("ABT644_Short_Gapping_List_Matches_The_Explicit_Spelling",
          "[constructive-model][core][gapping]") {
    // The one-entry form must be equivalent to spelling the residual gaps out by hand -- that
    // equivalence is the whole point, and it is what the bug broke. Checked on the derived gap
    // geometry, not just the lengths, and against the resulting magnetising inductance.
    auto build = [](json gapping) {
        json coreJson;
        coreJson["name"] = "gapping equivalence";
        coreJson["functionalDescription"] = json();
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "3C95";
        coreJson["functionalDescription"]["shape"] = "PQ 26/25";
        coreJson["functionalDescription"]["numberStacks"] = 1;
        coreJson["functionalDescription"]["gapping"] = gapping;
        return OpenMagnetics::Core(coreJson);
    };

    auto implicitCore = build(json::array({json{{"type", "subtractive"}, {"length", 0.0005}}}));
    auto explicitCore = build(json::array({
        json{{"type", "subtractive"}, {"length", 0.0005}},
        json{{"type", "residual"}, {"length", Constants().residualGap}},
        json{{"type", "residual"}, {"length", Constants().residualGap}}}));

    auto a = implicitCore.get_functional_description().get_gapping();
    auto b = explicitCore.get_functional_description().get_gapping();
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        UNSCOPED_INFO("gap index: " << i);
        CHECK(a[i].get_type() == b[i].get_type());
        CHECK_THAT(a[i].get_length(), Catch::Matchers::WithinRel(b[i].get_length(), 1e-9));
        CHECK_THAT(a[i].get_area().value(),
                   Catch::Matchers::WithinRel(b[i].get_area().value(), 1e-9));
        // Exact now. These two spellings land in different branches of
        // distribute_and_process_gap, and the branches used to disagree here: one subtracted half
        // the residual gap's length, the other ignored it (0.0080475 vs 0.0080500 on PQ 26/25).
        // Both now use (columnHeight - gapLength)/2, which is what the quantity means. ABT #644.
        CHECK_THAT(a[i].get_distance_closest_normal_surface().value(),
                   Catch::Matchers::WithinRel(b[i].get_distance_closest_normal_surface().value(),
                                              1e-9));
        // Including the AXIAL position. These two spellings used to reach different placement
        // paths and disagree by half a gap in y -- one centred the lone gap on the mating plane,
        // the other put it in one half. A ground gap is always machined into a single half, so
        // both now place it there and the whole coordinate triple must match. ABT #644.
        for (size_t axis : {0, 1, 2}) {
            UNSCOPED_INFO("axis: " << axis);
            CHECK_THAT(a[i].get_coordinates().value()[axis],
                       Catch::Matchers::WithinAbs(b[i].get_coordinates().value()[axis], 1e-12));
        }
    }

    // And the inductance must agree. Before the fix the implicit form read ~0.59x the explicit one
    // on this shape, which is the error that made this visible in the first place.
    auto coil = OpenMagneticsTesting::get_quick_coil(
        std::vector<int64_t>({10}), std::vector<int64_t>({1}), "PQ 26/25");
    MagnetizingInductance magnetizingInductance;
    auto implicitInductance = magnetizingInductance
        .calculate_inductance_from_number_turns_and_gapping(implicitCore, coil)
        .get_magnetizing_inductance().get_nominal().value();
    auto explicitInductance = magnetizingInductance
        .calculate_inductance_from_number_turns_and_gapping(explicitCore, coil)
        .get_magnetizing_inductance().get_nominal().value();
    UNSCOPED_INFO("implicit: " << implicitInductance << "  explicit: " << explicitInductance);
    CHECK_THAT(implicitInductance, Catch::Matchers::WithinRel(explicitInductance, 1e-9));
}

TEST_CASE("ABT644_Gap_Type_Decides_Where_The_Gap_Goes", "[constructive-model][core][gapping]") {
    // The gapping array is a vocabulary, and TYPE -- not position, not count -- decides which
    // column a gap lands on. The factory helpers define it:
    //     create_ground_gapping(L, n)         1 SUBTRACTIVE + (n-1) RESIDUAL
    //     create_distributed_gapping(L, N, n) N SUBTRACTIVE + (n-1) RESIDUAL
    //     create_spacer_gapping(L, n)         n ADDITIVE
    // So SUBTRACTIVE is always ground out of the CENTRAL column -- several of them are one
    // distributed gap spaced down that column, never one gap per leg -- while ADDITIVE is a shim
    // between the halves and therefore separates EVERY column.
    auto constants = Constants();
    const double L = 0.0005;

    auto build = [](json gapping) {
        json coreJson;
        coreJson["name"] = "gap vocabulary";
        coreJson["functionalDescription"] = json();
        coreJson["functionalDescription"]["type"] = "two-piece set";
        coreJson["functionalDescription"]["material"] = "3C95";
        coreJson["functionalDescription"]["shape"] = "PQ 26/25";
        coreJson["functionalDescription"]["numberStacks"] = 1;
        coreJson["functionalDescription"]["gapping"] = gapping;
        return OpenMagnetics::Core(coreJson);
    };
    auto sub = [&](double l) { return json{{"type", "subtractive"}, {"length", l}}; };
    auto add = [&](double l) { return json{{"type", "additive"}, {"length", l}}; };
    auto res = [&]() { return json{{"type", "residual"}, {"length", constants.residualGap}}; };

    auto countOfType = [](std::vector<CoreGap> const& g, GapType t) {
        size_t n = 0;
        for (auto const& x : g) if (x.get_type() == t) ++n;
        return n;
    };

    SECTION("three subtractive gaps are a DISTRIBUTED gap in the central column, not one per leg") {
        auto gapping = build(json::array({sub(L), sub(L), sub(L)}))
                           .get_functional_description().get_gapping();
        // three in the centre plus a residual on each of the two return columns
        REQUIRE(gapping.size() == 5);
        CHECK(countOfType(gapping, GapType::SUBTRACTIVE) == 3);
        CHECK(countOfType(gapping, GapType::RESIDUAL) == 2);
        // all three subtractive gaps share the central column's x, and are spread in y
        std::vector<double> heights;
        for (auto const& g : gapping) {
            if (g.get_type() != GapType::SUBTRACTIVE) continue;
            CHECK_THAT(g.get_coordinates().value()[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
            heights.push_back(g.get_coordinates().value()[1]);
        }
        REQUIRE(heights.size() == 3);
        std::sort(heights.begin(), heights.end());
        CHECK(heights[0] < heights[1]);
        CHECK(heights[1] < heights[2]);
        CHECK_THAT(heights[1], Catch::Matchers::WithinAbs(0.0, 1e-12));

        // and it must agree with the hand-built form, which is what callers use today
        auto explicitForm = build(json::array({sub(L), sub(L), sub(L), res(), res()}))
                                .get_functional_description().get_gapping();
        REQUIRE(explicitForm.size() == gapping.size());
        for (size_t i = 0; i < gapping.size(); ++i) {
            UNSCOPED_INFO("gap index: " << i);
            CHECK(gapping[i].get_type() == explicitForm[i].get_type());
            CHECK_THAT(gapping[i].get_length(),
                       Catch::Matchers::WithinRel(explicitForm[i].get_length(), 1e-9));
            for (size_t axis = 0; axis < 3; ++axis) {
                CHECK_THAT(gapping[i].get_coordinates().value()[axis],
                           Catch::Matchers::WithinAbs(
                               explicitForm[i].get_coordinates().value()[axis], 1e-12));
            }
        }
    }

    SECTION("a spacer separates EVERY column, however few additive gaps are given") {
        // one additive gap is still a shim between the halves: all three columns are pushed apart
        auto shortForm = build(json::array({add(L)}))
                             .get_functional_description().get_gapping();
        auto fullForm = build(json::array({add(L), add(L), add(L)}))
                            .get_functional_description().get_gapping();
        REQUIRE(shortForm.size() == 3);
        REQUIRE(fullForm.size() == 3);
        CHECK(countOfType(shortForm, GapType::ADDITIVE) == 3);
        CHECK(countOfType(fullForm, GapType::ADDITIVE) == 3);
        CHECK(countOfType(shortForm, GapType::RESIDUAL) == 0);
        for (size_t i = 0; i < shortForm.size(); ++i) {
            UNSCOPED_INFO("gap index: " << i);
            CHECK(shortForm[i].get_type() == fullForm[i].get_type());
            CHECK_THAT(shortForm[i].get_length(),
                       Catch::Matchers::WithinRel(fullForm[i].get_length(), 1e-9));
        }
    }

    SECTION("one subtractive gap is a ground gap: centre only, residual on the laterals") {
        auto gapping = build(json::array({sub(L)})).get_functional_description().get_gapping();
        REQUIRE(gapping.size() == 3);
        CHECK(countOfType(gapping, GapType::SUBTRACTIVE) == 1);
        CHECK(countOfType(gapping, GapType::RESIDUAL) == 2);
        CHECK(gapping[0].get_type() == GapType::SUBTRACTIVE);
        CHECK_THAT(gapping[0].get_coordinates().value()[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
        // ground into ONE half -- cheaper than machining both by half the gap each -- so the gap
        // spans 0..L and sits half a gap off the mating plane, never straddling it
        CHECK_THAT(gapping[0].get_coordinates().value()[1],
                   Catch::Matchers::WithinRel(L / 2, 1e-9));
        auto columnHeight = build(json::array({sub(L)}))
                                .get_processed_description().value().get_columns()[0].get_height();
        CHECK_THAT(gapping[0].get_distance_closest_normal_surface().value(),
                   Catch::Matchers::WithinRel(columnHeight / 2 - L / 2, 1e-6));
    }

    SECTION("an all-residual list stays one residual gap per column") {
        auto gapping = build(json::array({res(), res(), res()}))
                           .get_functional_description().get_gapping();
        REQUIRE(gapping.size() == 3);
        CHECK(countOfType(gapping, GapType::RESIDUAL) == 3);
    }
}

TEST_CASE("Toroid_Coating_Is_Resolved_Per_Material_Family", "[core][coating][abt964]") {
    // Which jacket a toroid carries, and how thick, is a property of the ceramic it is pressed
    // from. Epoxy on a FERRITE ring core is ~3x the film on a powder toroid, and the two families
    // cross over from parylene to epoxy at DIFFERENT sizes -- TDK puts a ferrite ring core's at
    // R 9.53 and Magnetics' ferrite catalogue agrees (coating code Y up to 7.62 mm OD, code Z
    // from 9.53 mm), while Micrometals/Fair-Rite put a powder toroid's at 0.20". Both thickness
    // numbers are tolerance-free: they compare a coated limit against the UNCOATED limit in the
    // same direction, so the bare core's dimensional tolerance cancels instead of being counted
    // as coating. TDK R 50.0x30.0x20.0 yields 0.400 mm on all three axes that way, matching the
    // "< 0.4 mm" it states; Ferroxcube draws its TN jacket at ~0.3 mm.
    auto toroid = [](const std::string& shape, const std::string& material) {
        json coreJson;
        coreJson["functionalDescription"]["type"] = "toroidal";
        coreJson["functionalDescription"]["material"] = material;
        coreJson["functionalDescription"]["shape"] = shape;
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        return Core(coreJson, true);
    };
    auto defaults = Defaults();

    // Ferrite past TDK's R 9.53 crossover: epoxy, at the drawn ring-core thickness.
    auto ferriteLarge = toroid("T 25/15/10", "N97");
    REQUIRE(ferriteLarge.is_ferrite_core());
    REQUIRE_FALSE(ferriteLarge.get_default_toroid_coating_is_parylene());
    REQUIRE_THAT(ferriteLarge.get_coating_thickness(),
                 Catch::Matchers::WithinAbs(defaults.defaultFerriteEpoxyCoreCoatingThickness, 1e-12));

    // THE CASE A SINGLE SHARED THRESHOLD GOT WRONG. One shape, 6,3 mm across, sits between the
    // two crossovers: past 0.20" so a POWDER toroid of that size is epoxy, but short of R 9.53 so
    // a FERRITE one is still parylene. Same geometry, opposite jackets, ~24x apart in thickness.
    auto ferriteMid = toroid("T 6.3/3.8/2.5", "N97");
    REQUIRE(ferriteMid.is_ferrite_core());
    REQUIRE(ferriteMid.get_default_toroid_coating_is_parylene());
    REQUIRE_THAT(ferriteMid.get_coating_thickness(),
                 Catch::Matchers::WithinAbs(defaults.defaultParyleneCoreCoatingThickness, 1e-12));

    auto powderMid = toroid("T 6.3/3.8/2.5", "Kool M\u00b5 26");
    REQUIRE_FALSE(powderMid.is_ferrite_core());
    REQUIRE_FALSE(powderMid.get_default_toroid_coating_is_parylene());
    REQUIRE_THAT(powderMid.get_coating_thickness(),
                 Catch::Matchers::WithinAbs(defaults.defaultEpoxyCoreCoatingThickness, 1e-12));

    REQUIRE(ferriteMid.get_coating_thickness() < powderMid.get_coating_thickness());
    REQUIRE(powderMid.get_coating_thickness() < ferriteLarge.get_coating_thickness());
}

TEST_CASE("Unprocessed toroid still resolves its ring edge", "[core][coating][abt964][regression]") {
    // REGRESSION. get_toroid_edge_radius() read the processed columns through
    //     const auto& columns = get_processed_description()->get_columns();
    // and get_processed_description() hands the optional back BY VALUE, so that reference pointed
    // into a temporary that was already gone. Undefined behaviour does not fail loudly: this build
    // read it fine, while asgard's WASM read it back EMPTY and threw "Toroid has no processed
    // column" for every catalogue toroid -- all 236 CMCs and El Choker's graph panel with them.
    //
    // The path that exposed it is the one autocomplete takes: a MAS record from the catalogue
    // carries no processedDescription, so the core arrives unprocessed and create_quick_bobbin
    // processes it on the way through. Build it that way here rather than pre-processing.
    json coreJson;
    coreJson["functionalDescription"]["type"] = "toroidal";
    coreJson["functionalDescription"]["material"] = "N97";
    coreJson["functionalDescription"]["shape"] = "T 25/15/10";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core unprocessed(coreJson, true);
    unprocessed.set_processed_description(std::nullopt);
    REQUIRE_FALSE(unprocessed.get_processed_description());

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(unprocessed, false);
    REQUIRE(bobbin.get_processed_description());
    REQUIRE(bobbin.get_column_corner_radius() > 0);

    // And the same core, processed up front, must agree: the datum is a property of the part,
    // not of the order the caller happened to build it in.
    Core processed(coreJson, true);
    auto processedBobbin = OpenMagnetics::Bobbin::create_quick_bobbin(processed, false);
    REQUIRE_THAT(bobbin.get_column_corner_radius(),
                 Catch::Matchers::WithinRel(processedBobbin.get_column_corner_radius(), 1e-12));
}
