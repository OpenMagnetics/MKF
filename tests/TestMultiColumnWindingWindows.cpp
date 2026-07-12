#include "constructive_models/Core.h"
#include "constructive_models/Bobbin.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

// Phase 1 of multi-column winding support: with the per-column-winding-windows setting
// on, Core::process_data emits one winding-window entry per wound-column candidate,
// each carrying the explicit window->column edge; with it off (default) the processed
// description is identical to the historical single-window output.

TEST_CASE("MultiColumnWindows_DefaultOff_HistoricalSingleWindow", "[constructive-model][core][multi-column][smoke-test]") {
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto windingWindows = core.get_processed_description()->get_winding_windows();
    REQUIRE(windingWindows.size() == 1);
    REQUIRE(!windingWindows[0].get_column());
    REQUIRE(core.get_main_column_index() == 0);
    // Schema-documented default: no column edge resolves to the main column.
    REQUIRE(core.get_winding_window_column_index(0) == 0);
    REQUIRE_THROWS(core.get_winding_window_column_index(1));
}

TEST_CASE("MultiColumnWindows_ECore_OneWindowPerColumn", "[constructive-model][core][multi-column][smoke-test]") {
    auto baselineCore = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto baselineWindow = baselineCore.get_processed_description()->get_winding_windows()[0];

    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto processedDescription = core.get_processed_description().value();
    auto windingWindows = processedDescription.get_winding_windows();
    auto columns = processedDescription.get_columns();

    REQUIRE(columns.size() == 3);
    REQUIRE(windingWindows.size() == 3);

    // Window 0 keeps its historical geometry and gains the main-column edge.
    CHECK(windingWindows[0].get_width().value() == baselineWindow.get_width().value());
    CHECK(windingWindows[0].get_height().value() == baselineWindow.get_height().value());
    CHECK(windingWindows[0].get_area().value() == baselineWindow.get_area().value());
    CHECK(windingWindows[0].get_coordinates().value()[0] == baselineWindow.get_coordinates().value()[0]);
    REQUIRE(windingWindows[0].get_column());
    CHECK(windingWindows[0].get_column().value() == 0);
    CHECK(columns[0].get_type() == ColumnType::CENTRAL);

    // One extra window per lateral column: same window geometry, placed on the
    // column's side of the main column, carrying that column's index.
    for (size_t windowIndex = 1; windowIndex < windingWindows.size(); ++windowIndex) {
        REQUIRE(windingWindows[windowIndex].get_column());
        auto columnIndex = static_cast<size_t>(windingWindows[windowIndex].get_column().value());
        REQUIRE(columnIndex < columns.size());
        CHECK(columns[columnIndex].get_type() == ColumnType::LATERAL);
        CHECK(windingWindows[windowIndex].get_width().value() == baselineWindow.get_width().value());
        CHECK(windingWindows[windowIndex].get_height().value() == baselineWindow.get_height().value());
        double columnX = columns[columnIndex].get_coordinates()[0];
        double windowX = windingWindows[windowIndex].get_coordinates().value()[0];
        CHECK(windowX * columnX > 0);
        CHECK(core.get_winding_window_column_index(windowIndex) == columnIndex);
    }

    // The two lateral windows must sit on opposite sides.
    CHECK(windingWindows[1].get_coordinates().value()[0] == -windingWindows[2].get_coordinates().value()[0]);
}

TEST_CASE("MultiColumnWindows_UCore_SharedRegionTwoColumns", "[constructive-model][core][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("U 15/11/6", json::parse("[]"), 1, "Dummy");
    auto processedDescription = core.get_processed_description().value();
    auto windingWindows = processedDescription.get_winding_windows();
    auto columns = processedDescription.get_columns();

    REQUIRE(columns.size() == 2);
    REQUIRE(windingWindows.size() == 2);
    REQUIRE(windingWindows[0].get_column());
    REQUIRE(windingWindows[1].get_column());
    CHECK(core.get_winding_window_column_index(0) == core.get_main_column_index());
    CHECK(static_cast<size_t>(windingWindows[1].get_column().value()) != core.get_main_column_index());
    // The U window region is shared by both legs: same geometry for both entries.
    CHECK(windingWindows[1].get_width().value() == windingWindows[0].get_width().value());
    CHECK(windingWindows[1].get_height().value() == windingWindows[0].get_height().value());
}

TEST_CASE("MultiColumnWindows_Toroid_SingleColumnStamped", "[constructive-model][core][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("T 17.3/9.7/12.7", json::parse("[]"), 1, "Dummy");
    auto windingWindows = core.get_processed_description()->get_winding_windows();
    REQUIRE(windingWindows.size() == 1);
    REQUIRE(windingWindows[0].get_column());
    CHECK(windingWindows[0].get_column().value() == 0);
    CHECK(core.get_winding_window_column_index(0) == 0);
}

TEST_CASE("MultiColumnWindows_QuickBobbin_PropagatesColumnEdge", "[constructive-model][core][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.0005, 0.0005);
    auto bobbinWindingWindows = bobbin.get_processed_description()->get_winding_windows();
    auto coreWindingWindows = core.get_processed_description()->get_winding_windows();

    REQUIRE(bobbinWindingWindows.size() == coreWindingWindows.size());
    for (size_t windowIndex = 0; windowIndex < bobbinWindingWindows.size(); ++windowIndex) {
        REQUIRE(bobbinWindingWindows[windowIndex].get_column());
        CHECK(bobbinWindingWindows[windowIndex].get_column().value() == coreWindingWindows[windowIndex].get_column().value());
        // Bobbin windows must sit on the same side as their core window.
        CHECK(bobbinWindingWindows[windowIndex].get_coordinates().value()[0] *
                  coreWindingWindows[windowIndex].get_coordinates().value()[0] > 0);
    }
}
