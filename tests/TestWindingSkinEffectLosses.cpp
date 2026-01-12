#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    double maximumError = 0.01;

    TEST_CASE("Test_Skin_Depth_Wire_Material_Data_20C", "[physical-model][skin-losses][smoke-test]") {
        auto materialData = find_wire_material_by_name("copper");
        double frequency = 123000;
        double temperature = 20;

        auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
        auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth(materialData, frequency, temperature);
        double expectedSkinDepth = 186.09e-6;
        REQUIRE_THAT(skinDepth, Catch::Matchers::WithinAbs(expectedSkinDepth, expectedSkinDepth * maximumError));
    }

    TEST_CASE("Test_Skin_Depth_Wire_Name_20C", "[physical-model][skin-losses][smoke-test]") {
        double frequency = 123000;
        double temperature = 20;

        auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
        auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);
        double expectedSkinDepth = 186.09e-6;
        REQUIRE_THAT(skinDepth, Catch::Matchers::WithinAbs(expectedSkinDepth, expectedSkinDepth * maximumError));
    }

    TEST_CASE("Test_Skin_Depth_Wire_Material_Data_200C", "[physical-model][skin-losses][smoke-test]") {
        auto materialData = find_wire_material_by_name("copper");
        double frequency = 123000;
        double temperature = 120;

        auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
        auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth(materialData, frequency, temperature);
        double expectedSkinDepth = 220e-6;
        REQUIRE_THAT(skinDepth, Catch::Matchers::WithinAbs(expectedSkinDepth, expectedSkinDepth * maximumError));
    }
// End of SUITE

}  // namespace
