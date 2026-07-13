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

TEST_CASE("MultiColumnWinding_ECore_PrimaryCentral_SecondaryLateral", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    auto columns = core.get_processed_description()->get_columns();
    auto coreWindingWindows = core.get_processed_description()->get_winding_windows();
    REQUIRE(coreWindingWindows.size() == 3);

    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 20}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
    coilJson["functionalDescription"].push_back(json{{"name", "Secondary"}, {"numberTurns", 10}, {"numberParallels", 1},
                                                     {"isolationSide", "secondary"}, {"wire", "Round 0.475 - Grade 1"}});

    OpenMagnetics::Coil coil(coilJson, false);
    // Place the secondary in the left lateral window (index 2: mirrored region, left lateral column).
    coil.get_mutable_functional_description()[1].set_winding_window(2);
    coil.set_core_columns(columns);
    REQUIRE(coil.wind());

    REQUIRE(coil.get_sections_description());
    REQUIRE(coil.get_layers_description());
    REQUIRE(coil.get_turns_description());

    size_t lateralColumnIndex = static_cast<size_t>(coreWindingWindows[2].get_column().value());
    auto lateralColumn = columns[lateralColumnIndex];
    REQUIRE(lateralColumn.get_type() == ColumnType::LATERAL);
    REQUIRE(lateralColumn.get_coordinates()[0] < 0);

    // Sections: primary in window 0 (positive x), secondary mirrored into the left window.
    auto woundSections = coil.get_sections_description().value();
    for (auto& section : woundSections) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        auto windingName = section.get_partial_windings()[0].get_winding();
        if (windingName == "Primary") {
            CHECK(section.get_coordinates()[0] > 0);
        }
        else {
            REQUIRE(section.get_winding_window());
            CHECK(section.get_winding_window().value() == 2);
            CHECK(section.get_coordinates()[0] < 0);
        }
    }

    // Turns: all present, secondary turns wrap the lateral column with the matching perimeter.
    size_t primaryTurns = 0;
    size_t secondaryTurns = 0;
    double lateralColumnAxisX = std::abs(lateralColumn.get_coordinates()[0]);
    // The lateral bobbin wraps its leg like the main one wraps the main column: one
    // column-wall thickness (0.001 in this quick bobbin) sits between leg and wire.
    double lateralBobbinColumnThickness = 0.001;
    double lateralHalfWidth = lateralColumn.get_width() / 2 + lateralBobbinColumnThickness;
    double lateralHalfDepth = lateralColumn.get_depth() / 2 + lateralBobbinColumnThickness;
    auto woundTurns = coil.get_turns_description().value();
    for (auto& turn : woundTurns) {
        if (turn.get_winding() == "Primary") {
            primaryTurns++;
            CHECK(turn.get_coordinates()[0] > 0);
            // Region sharing: the main winding's annulus keeps the inner half of the
            // window, so the lateral winding's outer half can never overlap it.
            CHECK(turn.get_coordinates()[0] < std::abs(coreWindingWindows[0].get_coordinates().value()[0]));
            CHECK(turn.get_length() > 0);
            // Second crossing of the main-column turn: the far side of the center leg.
            REQUIRE(turn.get_additional_coordinates());
            CHECK_THAT(turn.get_additional_coordinates().value()[0][0],
                       Catch::Matchers::WithinRel(-turn.get_coordinates()[0], 1e-9));
        }
        else {
            secondaryTurns++;
            CHECK(turn.get_coordinates()[0] < 0);
            // The winding hugs its leg: the turn sits on the leg side of the window.
            CHECK(std::abs(turn.get_coordinates()[0]) > std::abs(coreWindingWindows[2].get_coordinates().value()[0]));
            double radius = std::abs(std::abs(turn.get_coordinates()[0]) - lateralColumnAxisX);
            double expectedLength = 4 * lateralHalfDepth + 4 * lateralHalfWidth + 2 * std::numbers::pi * (radius - lateralHalfWidth);
            CHECK_THAT(turn.get_length(), Catch::Matchers::WithinRel(expectedLength, 1e-6));
            // Mirrored turns are wound the opposite way around their column.
            REQUIRE(turn.get_orientation());
            CHECK(turn.get_orientation().value() == TurnOrientation::COUNTER_CLOCKWISE);
            // Second crossing of the lateral turn: outside the core, mirrored across the leg axis.
            REQUIRE(turn.get_additional_coordinates());
            double secondCrossingX = turn.get_additional_coordinates().value()[0][0];
            CHECK(secondCrossingX < -core.get_processed_description()->get_width() / 2 * 0.99);
            CHECK_THAT(secondCrossingX, Catch::Matchers::WithinRel(2 * lateralColumn.get_coordinates()[0] - turn.get_coordinates()[0], 1e-9));
        }
    }
    CHECK(primaryTurns == 20);
    CHECK(secondaryTurns == 10);
}

TEST_CASE("MultiColumnWinding_DefaultPlacement_MatchesSingleWindow", "[constructive-model][coil][multi-column][smoke-test]") {
    // With the flag on but every winding left at the default placement (window 0),
    // the wound geometry must match the historical single-window result.
    auto buildAndWind = [](bool perColumnWindows) {
        auto& settings = Settings::GetInstance();
        settings.set_core_per_column_winding_windows(perColumnWindows);
        auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
        auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
        json coilJson;
        json bobbinJson;
        to_json(bobbinJson, bobbin);
        coilJson["bobbin"] = bobbinJson;
        coilJson["functionalDescription"] = json::array();
        coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 20}, {"numberParallels", 1},
                                                         {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
        coilJson["functionalDescription"].push_back(json{{"name", "Secondary"}, {"numberTurns", 10}, {"numberParallels", 1},
                                                         {"isolationSide", "secondary"}, {"wire", "Round 0.475 - Grade 1"}});
        OpenMagnetics::Coil coil(coilJson, false);
        coil.set_core_columns(core.get_processed_description()->get_columns());
        REQUIRE(coil.wind());
        return coil;
    };

    auto baselineCoil = buildAndWind(false);
    auto multiWindowCoil = buildAndWind(true);

    auto baselineTurns = baselineCoil.get_turns_description().value();
    auto multiWindowTurns = multiWindowCoil.get_turns_description().value();
    REQUIRE(baselineTurns.size() == multiWindowTurns.size());
    for (size_t turnIndex = 0; turnIndex < baselineTurns.size(); ++turnIndex) {
        CHECK(baselineTurns[turnIndex].get_coordinates()[0] == multiWindowTurns[turnIndex].get_coordinates()[0]);
        CHECK(baselineTurns[turnIndex].get_coordinates()[1] == multiWindowTurns[turnIndex].get_coordinates()[1]);
        CHECK(baselineTurns[turnIndex].get_length() == multiWindowTurns[turnIndex].get_length());
    }
}

