#include <source_location>
#include "constructive_models/Mas.h"
#include "support/Settings.h"
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
    TEST_CASE("Test_Expand_Magnetic", "[constructive-model][mas][smoke-test]") {

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "example_basic.json");
        auto mas = OpenMagneticsTesting::mas_loader(path.string());

        REQUIRE(!mas.get_magnetic().get_core().get_processed_description());

        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
        REQUIRE(magnetic.get_core().get_processed_description());

    }

}  // namespace
