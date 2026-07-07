#include <source_location>
#include "constructive_models/Bobbin.h"
#include "constructive_models/Core.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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
const auto masPath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("MAS/").string();
const double max_error = 0.05;

// ABT #107: every named-catalogue bobbin family must store the winding-window
// x-coordinate as the window CENTRE (consumers compute left = coords[0] - width/2,
// matching create_quick_bobbin). This regression test locks that convention across
// both historical bug classes so the three coexisting conventions can never return:
//   - inner-edge families (E/EP/PM/RM/EFD) that stored only the column surface, and
//   - full-diameter families (ETD/PQ/EC) that stored 2× the inner radius (the centre
//     ended up beyond the window's own outer edge).
// The invariant checked for each: coords[0] - width/2 == the column-surface inner
// edge (positive, sane), i.e. coords[0] is the true centre. Values captured from the
// datasheet dimensions of the named bobbins below.
TEST_CASE("Bobbin winding window coordinate is the window centre (ABT #107)",
          "[constructive-model][bobbin][abt107]") {
    struct Case { std::string name; double expectedInnerEdge; };
    // RM 10  (round, inner-edge family): D2/2 = 0.0062, width (D1-D2)/2 = 0.00425
    // EC 35  (round, full-diameter family): D2/2 = 0.0061, width (D1-D2)/2 = 0.004675
    const std::vector<Case> cases = {
        {"Bobbin RM 10", 0.0062},
        {"Bobbin EC 35", 0.0061},
    };
    for (const auto& c : cases) {
        auto bobbin = OpenMagnetics::find_bobbin_by_name(c.name);
        auto proc = bobbin.get_processed_description().value();
        auto windingWindows = proc.get_winding_windows();
        REQUIRE(!windingWindows.empty());
        REQUIRE(windingWindows[0].get_coordinates().has_value());
        REQUIRE(windingWindows[0].get_width().has_value());
        double centre = windingWindows[0].get_coordinates().value()[0];
        double width = windingWindows[0].get_width().value();
        double leftEdge = centre - width / 2;
        INFO("bobbin=" << c.name << " centre=" << centre << " width=" << width << " leftEdge=" << leftEdge);
        // Left edge sits at the column surface, not inside the column or at 2× it.
        REQUIRE_THAT(leftEdge, Catch::Matchers::WithinRel(c.expectedInnerEdge, 1e-6));
        // Centre is strictly outside the column surface (the historical inner-edge bug
        // put the centre AT the surface; the full-diameter bug put it past the far edge).
        REQUIRE(centre > leftEdge);
        REQUIRE(centre < leftEdge + width);
    }
}

// The E19_5 sample (rectangular E family) exercises the third processor path.
TEST_CASE("Bobbin E sample winding window is centred (ABT #107)",
          "[constructive-model][bobbin][abt107]") {
    auto wireFilePath = masPath + "samples/magnetic/bobbin/bobbin_E19_5.json";
    std::ifstream json_file(wireFilePath);
    auto bobbinJson = json::parse(json_file);
    OpenMagnetics::Bobbin bobbin(bobbinJson);
    auto proc = bobbin.get_processed_description().value();
    auto windingWindows = proc.get_winding_windows();
    double centre = windingWindows[0].get_coordinates().value()[0];
    double width = windingWindows[0].get_width().value();
    // Inner edge = f/2 + s1 = 0.00485/2 + 0.0008 = 0.003225 (column surface).
    REQUIRE_THAT(centre - width / 2, Catch::Matchers::WithinRel(0.003225, 1e-6));
}

// ABT #107: the CORE winding-window coordinate must also be the window CENTRE.
// The E-family core pieces stored only the inner edge (F/2); consumers (PainterImpl,
// create_quick_bobbin sign check) and the U/Ur/C pieces expect the centre. Invariant:
// window left edge (coords[0] - width/2) == central-column half-width (the column surface).
TEST_CASE("Core winding window coordinate is the window centre (ABT #107)",
          "[constructive-model][core][abt107]") {
    auto core = OpenMagnetics::Core::create_quick_core("E 42/21/15", "3C97");
    auto processed = core.get_processed_description().value();
    auto windingWindows = processed.get_winding_windows();
    auto columns = processed.get_columns();
    REQUIRE(!windingWindows.empty());
    REQUIRE(!columns.empty());
    double centre = windingWindows[0].get_coordinates().value()[0];
    double width = windingWindows[0].get_width().value();
    double columnHalfWidth = columns[0].get_width() / 2;   // F/2, the central-column surface
    INFO("centre=" << centre << " width=" << width << " columnHalfWidth=" << columnHalfWidth);
    REQUIRE_THAT(centre - width / 2, Catch::Matchers::WithinRel(columnHalfWidth, 1e-6));
    REQUIRE(centre > columnHalfWidth);   // strictly outside the column surface (was AT it before the fix)
}

TEST_CASE("Sample_Bobbin", "[constructive-model][bobbin][smoke-test]") {
    auto wireFilePath = masPath + "samples/magnetic/bobbin/bobbin_E19_5.json";
    std::ifstream json_file(wireFilePath);
    auto bobbinJson = json::parse(json_file);

    OpenMagnetics::Bobbin bobbin(bobbinJson);

    double expectedColumnThickness = 0.00080;
    double expectedWallThickness = 0.00080;

    REQUIRE_THAT(expectedColumnThickness, Catch::Matchers::WithinAbs(bobbin.get_processed_description().value().get_column_thickness(), max_error * expectedColumnThickness));
    REQUIRE_THAT(expectedWallThickness, Catch::Matchers::WithinAbs(bobbin.get_processed_description().value().get_wall_thickness(), max_error * expectedWallThickness));
}

TEST_CASE("Get_Filling_Factors_Bobbin_Medium", "[constructive-model][bobbin][smoke-test]") {
    auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.009, 0.0275);

    double expectedValue = 0.715;

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
}

TEST_CASE("Get_Filling_Factors_Bobbin_Small", "[constructive-model][bobbin][smoke-test]") {
    auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.002, 0.005);

    double expectedValue = 0.53;

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
}

TEST_CASE("Get_Filling_Factors_Bobbin_Large", "[constructive-model][bobbin][smoke-test]") {
    auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.019, 0.057);

    double expectedValue = 0.725;

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
}

TEST_CASE("Get_Filling_Factors_Bobbin_Outside_Above", "[constructive-model][bobbin][smoke-test]") {
    auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(1, 1);

    double expectedValue = 0.79;

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
}

TEST_CASE("Get_Filling_Factors_Bobbin_Outside_Below", "[constructive-model][bobbin][smoke-test]") {
    auto fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0, 0);

    double expectedValue = 0.377;

    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
}

TEST_CASE("Get_Winding_Window_Dimensions_Medium", "[constructive-model][bobbin][smoke-test]") {
    auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.012, 0.027);

    double expectedWidthValue = 0.00985;
    double expectedHeightValue = 0.02335;
    double width = windingWindowDimensions[0];
    double height = windingWindowDimensions[1];

    REQUIRE_THAT(expectedWidthValue, Catch::Matchers::WithinAbs(width, max_error * expectedWidthValue));
    REQUIRE_THAT(expectedHeightValue, Catch::Matchers::WithinAbs(height, max_error * expectedHeightValue));
}

TEST_CASE("Get_Winding_Window_Dimensions_Too_Small", "[constructive-model][bobbin][smoke-test]") {
    auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.001, 0.002);

    double expectedWidthValue = 0.0005;
    double expectedHeightValue = 0.0011;
    double width = windingWindowDimensions[0];
    double height = windingWindowDimensions[1];

    REQUIRE_THAT(expectedWidthValue, Catch::Matchers::WithinAbs(width, max_error * expectedWidthValue));
    REQUIRE_THAT(expectedHeightValue, Catch::Matchers::WithinAbs(height, max_error * expectedHeightValue));
}

TEST_CASE("Get_Winding_Window_Dimensions_Too_Large", "[constructive-model][bobbin][smoke-test]") {
    auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.1, 0.1);

    double expectedWidthValue = 0.0951;
    double expectedHeightValue = 0.0943;
    double width = windingWindowDimensions[0];
    double height = windingWindowDimensions[1];

    REQUIRE_THAT(expectedWidthValue, Catch::Matchers::WithinAbs(width, max_error * expectedWidthValue));
    REQUIRE_THAT(expectedHeightValue, Catch::Matchers::WithinAbs(height, max_error * expectedHeightValue));
}

TEST_CASE("Get_Winding_Window_Dimensions_Error", "[constructive-model][bobbin][smoke-test]") {
    auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(0.003325, 0.0108);

    double expectedWidthValue =  0.00245;
    double expectedHeightValue = 0.0094;
    double width = windingWindowDimensions[0];
    double height = windingWindowDimensions[1];

    REQUIRE_THAT(expectedWidthValue, Catch::Matchers::WithinAbs(width, max_error * expectedWidthValue));
    REQUIRE_THAT(expectedHeightValue, Catch::Matchers::WithinAbs(height, max_error * expectedHeightValue));
}

TEST_CASE("Get_Winding_Window_Dimensions_E_51", "[constructive-model][bobbin][smoke-test]") {
    auto core = OpenMagneticsTesting::get_quick_core("ER 51/10/38", json::parse("[]"), 1, "Dummy");
    auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
    auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(coreWindingWindow.get_width().value(), coreWindingWindow.get_height().value());

    auto widthThickness = coreWindingWindow.get_width().value() - windingWindowDimensions[0];
    auto heightThickness = (coreWindingWindow.get_height().value() - windingWindowDimensions[1]) / 2;

    REQUIRE(heightThickness <= widthThickness * 1.2);
}

TEST_CASE("Get_Winding_Window_Dimensions_All_Shapes_With_Bobbin", "[constructive-model][bobbin][smoke-test]") {
    settings.set_use_toroidal_cores(true);
    auto shapeNames = get_core_shape_names();
    for (auto shapeName : shapeNames) {
        if (shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI ")) {
            continue;
        }
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
        auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
        auto windingWindowDimensions = OpenMagnetics::Bobbin::get_winding_window_dimensions(coreWindingWindow.get_width().value(), coreWindingWindow.get_height().value());

        // The bobbin winding window must be positive and fit inside the core winding window.
        INFO("Shape: " << shapeName);
        REQUIRE(windingWindowDimensions.size() == 2);
        CHECK(windingWindowDimensions[0] > 0);
        CHECK(windingWindowDimensions[0] <= coreWindingWindow.get_width().value());
        CHECK(windingWindowDimensions[1] > 0);
        CHECK(windingWindowDimensions[1] <= coreWindingWindow.get_height().value());
    }
}

TEST_CASE("Create_Bobbin_With_Thickness", "[constructive-model][bobbin][smoke-test]") {
    settings.set_use_toroidal_cores(true);
    auto shapeNames = get_core_shape_names();
    for (auto shapeName : shapeNames) {
        if (shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI ")) {
            continue;
        }
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
        auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
        double wallThickness = coreWindingWindow.get_height().value() * 0.1;
        double columnThickness = coreWindingWindow.get_width().value() * 0.1;
        auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, wallThickness, columnThickness);
        REQUIRE_THAT(wallThickness, Catch::Matchers::WithinAbs(bobbin.get_processed_description().value().get_wall_thickness(), max_error * wallThickness));
        REQUIRE_THAT(columnThickness, Catch::Matchers::WithinAbs(bobbin.get_processed_description().value().get_column_thickness(), max_error * columnThickness));
    }
}

}  // namespace
