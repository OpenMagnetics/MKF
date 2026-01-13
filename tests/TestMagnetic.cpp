#include "physical_models/MagneticEnergy.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

using json = nlohmann::json;

namespace {
    TEST_CASE("Test_Magnetic_Saturation_Current", "[constructive-model][magnetic][smoke-test]") {

        std::vector<int64_t> numberTurns = {18};
        std::vector<int64_t> numberParallels = {1};
        std::string shapeName = "PQ 65/44";

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName);

        // coil.wind({0, 1, 0}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto saturationCurrentAt20 = magnetic.calculate_saturation_current(20);
        auto saturationCurrentAt100 = magnetic.calculate_saturation_current(100);
        REQUIRE(saturationCurrentAt100 < saturationCurrentAt20);
    }

}  // namespace
