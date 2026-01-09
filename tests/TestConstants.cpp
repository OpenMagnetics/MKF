#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
using json = nlohmann::json;

namespace {

TEST_CASE("AccessConstants", "[support][constants]") {
    REQUIRE_THAT(OpenMagnetics::constants.residualGap, Catch::Matchers::WithinAbs(5.0e-6, 1e-8));
}

}  // namespace
