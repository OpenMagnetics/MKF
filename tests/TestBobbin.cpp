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
#include <thread>
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

// ABT #631: a bobbin catalogue row that points at a core shape MAS does not ship (34 of
// the 504 shipped bobbins do: 10 "EI …", 24 "M …") must be SKIPPED by the interpolator
// fit, and the skip must not be spelled as throw-and-catch. Wherever MKF is compiled with
// exception catching disabled — the Emscripten default, and how MVB++'s WASM module builds
// it — a catch handler is deleted outright, so a throw-and-catch skip becomes a fatal error
// that escapes to the caller: in the browser the first such row ("Bobbin EI 101/50") took
// down every design whose bobbin or core shape was given by NAME. The load must also SAY
// what it dropped; a bare `continue` is how 34 dangling references went unnoticed.
TEST_CASE("Unresolvable bobbin catalogue rows are skipped without throwing (ABT #631)",
          "[constructive-model][bobbin][abt631]") {
    SECTION("core shape lookup can be asked instead of thrown at") {
        CHECK(OpenMagnetics::core_shape_exists("ETD 49/25/16"));
        CHECK(OpenMagnetics::try_find_core_shape_by_name("ETD 49/25/16").has_value());
        // The name the browser build died on. It is a real bobbin's shape reference and
        // NOT a shape MAS ships; asking must answer "no", not throw.
        CHECK_FALSE(OpenMagnetics::core_shape_exists("EI 101/50"));
        CHECK_FALSE(OpenMagnetics::try_find_core_shape_by_name("EI 101/50").has_value());
        CHECK_THROWS_AS(OpenMagnetics::find_core_shape_by_name("EI 101/50"), OpenMagnetics::CoreException);
    }

    SECTION("the shipped catalogue still carries the dangling references it always did") {
        auto bobbin = OpenMagnetics::find_bobbin_by_name("Bobbin EI 101/50");
        REQUIRE(bobbin.get_functional_description());
        CHECK(bobbin.get_functional_description()->get_shape() == "EI 101/50");
        CHECK_FALSE(OpenMagnetics::core_shape_exists(bobbin.get_functional_description()->get_shape()));
    }

    SECTION("an injected dangling row is skipped, reported, and does not break the fit") {
        // Poison the shared catalogue with a row pointing at a shape that cannot exist,
        // then fit the interpolators from scratch. The interpolators are thread_local
        // (ABT #113), so a FRESH thread is what guarantees a real fit here: this thread's
        // splines are already warm from the tests above and would not be refitted.
        auto poison = OpenMagnetics::find_bobbin_by_name("Bobbin ETD 49");
        auto functionalDescription = poison.get_functional_description().value();
        functionalDescription.set_shape("Shape That Does Not Exist 631");
        poison.set_functional_description(functionalDescription);
        const std::string poisonName = "Bobbin ABT631 Poison";
        OpenMagnetics::bobbinDatabase[poisonName] = poison;
        OpenMagnetics::read_log();  // drain, so the assertions below see only this fit

        double fillingFactor = 0;
        std::string threadError;
        std::thread fitter([&] {
            try {
                fillingFactor = OpenMagnetics::Bobbin::get_filling_factor(0.01, 0.01);
            }
            catch (const std::exception& e) {
                threadError = e.what();
            }
        });
        fitter.join();
        OpenMagnetics::bobbinDatabase.erase(poisonName);

        // The whole point: one unresolvable row must not take the fit down with it.
        INFO("fit threw: " << threadError);
        CHECK(threadError.empty());
        CHECK(fillingFactor > 0);

        auto log = OpenMagnetics::read_log();
        // The skip is visible and names the row and the reason, so a dangling reference
        // can be found and fixed instead of silently shrinking the fit sample.
        CHECK(log.find(poisonName) != std::string::npos);
        CHECK(log.find("Shape That Does Not Exist 631") != std::string::npos);
    }

    SECTION("an unresolved placeholder bobbin says so instead of reading uninitialised memory") {
        // find_bobbin_by_name returns a bobbin with NEITHER description for the documented
        // placeholder names, and Coil::resolve_bobbin hands that straight to callers. This
        // used to dereference the disengaged functionalDescription optional — undefined
        // behaviour that read whatever the returned-by-value optional's storage held.
        for (const auto& placeholder : {"basic", "Basic", "Dummy", "None"}) {
            auto bobbin = OpenMagnetics::find_bobbin_by_name(placeholder);
            REQUIRE_FALSE(bobbin.get_processed_description());
            REQUIRE_FALSE(bobbin.get_functional_description());
            CHECK_THROWS_AS(bobbin.get_winding_window_shape(), OpenMagnetics::InvalidInputException);
        }
    }
}

}  // namespace
