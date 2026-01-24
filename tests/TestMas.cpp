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

    TEST_CASE("Test_Load_Mas_With_Microsign_In_Wire_Name", "[constructive-model][mas][encoding][smoke-test]") {
        // This test verifies that files containing the micro sign (µ) character in wire names
        // are loaded correctly across platforms (Windows/Linux encoding differences)
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_microsign_encoding.json");
        
        std::ifstream json_file(path);
        REQUIRE(json_file.is_open());
        
        auto jsonData = json::parse(json_file);
        json_file.close();
        
        // Verify the JSON was parsed successfully
        REQUIRE(jsonData.contains("magnetic"));
        REQUIRE(jsonData["magnetic"].contains("coil"));
        REQUIRE(jsonData["magnetic"]["coil"].contains("functionalDescription"));
        
        // Get the wire name from the first winding
        auto& windings = jsonData["magnetic"]["coil"]["functionalDescription"];
        REQUIRE(windings.size() > 0);
        REQUIRE(windings[0].contains("wire"));
        REQUIRE(windings[0]["wire"].contains("name"));
        
        std::string wireName = windings[0]["wire"]["name"];
        
        // The wire name should contain "Planar" and should be non-empty
        REQUIRE(!wireName.empty());
        REQUIRE(wireName.find("Planar") != std::string::npos);
        
        // Verify we can construct a Mas object from the JSON
        OpenMagnetics::Mas mas(jsonData);
        
        // Get the wire name from the loaded Mas object
        // get_wire() returns std::variant<Wire, std::string>
        const auto& wireVariant = mas.get_magnetic().get_coil().get_functional_description()[0].get_wire();
        REQUIRE(std::holds_alternative<OpenMagnetics::Wire>(wireVariant));
        const auto& wire = std::get<OpenMagnetics::Wire>(wireVariant);
        auto loadedWireName = wire.get_name();
        REQUIRE(loadedWireName.has_value());
        REQUIRE(loadedWireName.value().find("Planar") != std::string::npos);
        
        // Verify the magnetic can be processed
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        REQUIRE(magnetic.get_core().get_processed_description());
    }

}  // namespace
