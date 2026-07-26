#include "constructive_models/Core.h"
#include "constructive_models/Bobbin.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <fstream>
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

// ---------------------------------------------------------------------------
// JSON-boundary + custom-rectangle regressions (July 2026 winding-studio audit).
// Every bug below shipped silently: the C++ tests all built their coils in
// process, so nothing exercised the serializers or the hand-edited-section
// paths the web studio drives.
// ---------------------------------------------------------------------------

namespace {
struct WoundLateralE42 {
    Core core;
    std::vector<ColumnElement> columns;
    OpenMagnetics::Coil coil;
};

WoundLateralE42 build_wound_e42_lateral_secondary() {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    auto columns = core.get_processed_description()->get_columns();
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
    coil.get_mutable_functional_description()[1].set_winding_window(2);
    coil.set_core_columns(columns);
    REQUIRE(coil.wind());
    return {core, columns, coil};
}

double section_x(const OpenMagnetics::Coil& coil, const std::string& windingName) {
    auto sections = coil.get_sections_description().value();
    for (auto& section : sections) {
        if (section.get_type() == ElectricalType::CONDUCTION &&
            section.get_partial_windings()[0].get_winding() == windingName) {
            return section.get_coordinates()[0];
        }
    }
    throw std::runtime_error("no conduction section for " + windingName);
}
}  // namespace

// f8ff5789: OpenMagnetics::Winding's hand-rolled from_json/to_json dropped the
// winding-level windingWindow, so any placement set through a JSON boundary
// (WASM bindings, file loads) was silently wound into window 0.
TEST_CASE("MultiColumnWinding_JsonBoundary_WindingWindowSurvives", "[constructive-model][coil][multi-column][smoke-test]") {
    auto wound = build_wound_e42_lateral_secondary();

    json serialized;
    to_json(serialized, wound.coil);
    REQUIRE(serialized["functionalDescription"][1].contains("windingWindow"));
    CHECK(serialized["functionalDescription"][1]["windingWindow"] == 2);

    OpenMagnetics::Coil reloaded(serialized, false);
    REQUIRE(reloaded.get_functional_description()[1].get_winding_window());
    CHECK(reloaded.get_functional_description()[1].get_winding_window().value() == 2);

    // A fresh wind of the reloaded coil still places the secondary laterally.
    reloaded.set_core_columns(wound.columns);
    REQUIRE(reloaded.wind());
    CHECK(section_x(reloaded, "Primary") > 0);
    CHECK(section_x(reloaded, "Secondary") < 0);
}

// 5eb32261: a coil reconstructed from JSON starts with the placement-transform
// flag unset while its serialized sections are at their FINAL mirrored
// positions; re-compacting moved the lateral section back to the main window.
TEST_CASE("MultiColumnWinding_JsonReload_DelimitKeepsMirroredSections", "[constructive-model][coil][multi-column][smoke-test]") {
    auto wound = build_wound_e42_lateral_secondary();

    json serialized;
    to_json(serialized, wound.coil);
    OpenMagnetics::Coil reloaded(serialized, false);
    reloaded.set_core_columns(wound.columns);
    REQUIRE(section_x(reloaded, "Secondary") < 0);

    REQUIRE(reloaded.delimit_and_compact());
    CHECK(section_x(reloaded, "Secondary") < 0);
    auto reloadedTurns = reloaded.get_turns_description().value();
    for (auto& turn : reloadedTurns) {
        if (turn.get_winding() == "Secondary") {
            CHECK(turn.get_coordinates()[0] < 0);
        }
    }
}

// bdd2ae16 + e5f5d609 + 48d05d4c: a hand-drawn rectangle on a MIRRORED-window
// section re-flows layers+turns INSIDE the rectangle (through the frame
// round-trip), re-packing them (the stale numberLayers is cleared) and
// surviving subsequent full winds.
TEST_CASE("MultiColumnWinding_CustomSectionRect_RepacksAndSurvivesRewind", "[constructive-model][coil][multi-column][smoke-test]") {
    auto wound = build_wound_e42_lateral_secondary();
    auto& coil = wound.coil;

    auto sections = coil.get_sections_description().value();
    std::string sectionName;
    std::vector<double> baseCoordinates;
    std::vector<double> baseDimensions;
    for (auto& section : sections) {
        if (section.get_type() == ElectricalType::CONDUCTION &&
            section.get_partial_windings()[0].get_winding() == "Secondary") {
            sectionName = section.get_name();
            baseCoordinates = section.get_coordinates();
            baseDimensions = section.get_dimensions();
        }
    }
    REQUIRE(!sectionName.empty());

    // 2.2x wider (room for a second radial layer, clear of the exact-multiple
    // float truncation) and 55% of the height, pushed toward the window top.
    std::vector<double> customDimensions = {baseDimensions[0] * 2.2, baseDimensions[1] * 0.55};
    std::vector<double> customCoordinates = {baseCoordinates[0], 0.0283 / 2 - customDimensions[1] / 2};
    coil.preload_custom_section_rects({{sectionName, {customCoordinates, customDimensions}}});
    REQUIRE(coil.apply_custom_section_rects());

    auto check_custom_layout = [&]() {
        auto woundSections = coil.get_sections_description().value();
        size_t conductionLayers = 0;
        for (auto& section : woundSections) {
            if (section.get_name() == sectionName) {
                CHECK_THAT(section.get_coordinates()[0], Catch::Matchers::WithinRel(customCoordinates[0], 1e-9));
                CHECK_THAT(section.get_coordinates()[1], Catch::Matchers::WithinRel(customCoordinates[1], 1e-9));
                CHECK_THAT(section.get_dimensions()[1], Catch::Matchers::WithinRel(customDimensions[1], 1e-9));
            }
        }
        auto woundLayers = coil.get_layers_description().value();
        for (auto& layer : woundLayers) {
            if (layer.get_type() == ElectricalType::CONDUCTION && layer.get_section() == sectionName) {
                conductionLayers++;
            }
        }
        CHECK(conductionLayers == 2);
        double yLow = customCoordinates[1] - customDimensions[1] / 2;
        double yHigh = customCoordinates[1] + customDimensions[1] / 2;
        auto woundTurns = coil.get_turns_description().value();
        for (auto& turn : woundTurns) {
            if (turn.get_winding() == "Secondary") {
                CHECK(turn.get_coordinates()[0] < 0);
                CHECK(turn.get_coordinates()[1] > yLow - 1e-9);
                CHECK(turn.get_coordinates()[1] < yHigh + 1e-9);
            }
        }
    };
    check_custom_layout();

    // The pin survives a FULL rewind with different proportions (compaction
    // runs for the rest of the coil; the drawn section is re-imposed after).
    REQUIRE(coil.wind(std::vector<double>{0.7, 0.3}, std::vector<size_t>{0, 1}, 1));
    check_custom_layout();
}

// Per-column bobbin ARRAY (MAS coil.bobbin array form): a catalog-style centre
// bobbin plus an ad-hoc lateral bobbin — two BOM items. The boundary must
// merge them into one effective multi-window bobbin for winding while
// serialization keeps the honest two-element array. Each element's windows
// carry their own column edge, which is what the winder follows.
TEST_CASE("MultiColumnWinding_PerColumnBobbinArray_MergesWindsAndRoundTrips", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto quickBobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    auto columns = core.get_processed_description()->get_columns();
    auto quickWindows = quickBobbin.get_processed_description()->get_winding_windows();
    REQUIRE(quickWindows.size() == 3);

    // Centre part: window 0 only (what a real catalog bobbin describes).
    auto centerBobbin = quickBobbin;
    auto centerProcessedDescription = centerBobbin.get_processed_description().value();
    centerProcessedDescription.set_winding_windows({quickWindows[0]});
    centerBobbin.set_processed_description(centerProcessedDescription);

    // Ad-hoc lateral part: the left-lateral window only, its column edge intact.
    auto lateralBobbin = quickBobbin;
    auto lateralProcessedDescription = lateralBobbin.get_processed_description().value();
    lateralProcessedDescription.set_winding_windows({quickWindows[2]});
    lateralBobbin.set_processed_description(lateralProcessedDescription);
    size_t lateralColumnIndex = static_cast<size_t>(quickWindows[2].get_column().value());
    REQUIRE(columns[lateralColumnIndex].get_type() == ColumnType::LATERAL);

    json centerJson;
    json lateralJson;
    to_json(centerJson, centerBobbin);
    to_json(lateralJson, lateralBobbin);
    json coilJson;
    coilJson["bobbin"] = json::array({centerJson, lateralJson});
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 20}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
    coilJson["functionalDescription"].push_back(json{{"name", "Secondary"}, {"numberTurns", 10}, {"numberParallels", 1},
                                                     {"isolationSide", "secondary"}, {"wire", "Round 0.475 - Grade 1"}});

    OpenMagnetics::Coil coil(coilJson, false);

    // The effective bobbin is the merge: centre window + lateral window.
    auto mergedWindows = coil.resolve_bobbin().get_processed_description()->get_winding_windows();
    REQUIRE(mergedWindows.size() == 2);
    REQUIRE(mergedWindows[1].get_column());
    CHECK(static_cast<size_t>(mergedWindows[1].get_column().value()) == lateralColumnIndex);

    // Wind the secondary in the lateral window (merged index 1).
    coil.get_mutable_functional_description()[1].set_winding_window(1);
    coil.set_core_columns(columns);
    REQUIRE(coil.wind());
    CHECK(section_x(coil, "Primary") > 0);
    CHECK(section_x(coil, "Secondary") < 0);
    size_t secondaryTurns = 0;
    auto woundTurns = coil.get_turns_description().value();
    for (auto& turn : woundTurns) {
        if (turn.get_winding() == "Secondary") {
            secondaryTurns++;
            CHECK(turn.get_coordinates()[0] < 0);
        }
    }
    CHECK(secondaryTurns == 10);

    // Serialization must re-emit the two-element array (the honest BOM), with
    // each part keeping its own single window — never the merged working bobbin.
    json serialized;
    to_json(serialized, coil);
    REQUIRE(serialized["bobbin"].is_array());
    REQUIRE(serialized["bobbin"].size() == 2);
    CHECK(serialized["bobbin"][0]["processedDescription"]["windingWindows"].size() == 1);
    CHECK(serialized["bobbin"][1]["processedDescription"]["windingWindows"].size() == 1);

    // Reload and rewind: placement and merge survive the round trip.
    OpenMagnetics::Coil reloaded(serialized, false);
    REQUIRE(reloaded.get_per_column_bobbins());
    CHECK(reloaded.get_per_column_bobbins()->size() == 2);
    reloaded.set_core_columns(columns);
    CHECK(section_x(reloaded, "Secondary") < 0);
    REQUIRE(reloaded.wind());
    CHECK(section_x(reloaded, "Secondary") < 0);
}

// Winding-style override (winding studio): the user can force multifilar
// (consecutive parallels) or turn-by-turn (consecutive turns) per winding;
// the heuristic keeps deciding for windings without an override.
TEST_CASE("WindingStyleOverride_ForcesConsecutiveTurns", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(false);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 20}, {"numberParallels", 2},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});

    auto windWithOverride = [&](std::optional<WindingStyle> styleOverride) {
        OpenMagnetics::Coil coil(coilJson, false);
        if (styleOverride) {
            coil.preload_winding_style_overrides({{"Primary", styleOverride.value()}});
        }
        REQUIRE(coil.wind());
        auto sections = coil.get_sections_description().value();
        for (auto& section : sections) {
            if (section.get_type() == ElectricalType::CONDUCTION) {
                REQUIRE(section.get_winding_style());
                return section.get_winding_style().value();
            }
        }
        throw std::runtime_error("no conduction section");
    };

    // Heuristic default for 20 turns × 2 parallels in one slot: parallels grouped.
    CHECK(windWithOverride(std::nullopt) == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
    // Override flips it; the wind still lands every turn.
    CHECK(windWithOverride(WindingStyle::WIND_BY_CONSECUTIVE_TURNS) == WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
}

// woundWith through the JSON boundary: same silent-drop bug class as the
// winding-level windingWindow — the custom Winding serializer must carry it,
// or WASM/file loads never group the windings into shared sections.
TEST_CASE("WoundWith_SurvivesJsonBoundary_AndGroupsSections", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(false);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 12}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"},
                                                     {"woundWith", json::array({"Bias"})}});
    coilJson["functionalDescription"].push_back(json{{"name", "Bias"}, {"numberTurns", 12}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"},
                                                     {"woundWith", json::array({"Primary"})}});

    OpenMagnetics::Coil coil(coilJson, false);
    REQUIRE(coil.get_functional_description()[0].get_wound_with());
    CHECK(coil.get_functional_description()[0].get_wound_with().value() == std::vector<std::string>{"Bias"});
    REQUIRE(coil.wind());

    // The pair shares its conduction sections (bifilar: both partial windings
    // in the same section) and every turn of both windings lands.
    auto sections = coil.get_sections_description().value();
    bool sharedSectionFound = false;
    for (auto& section : sections) {
        if (section.get_type() == ElectricalType::CONDUCTION && section.get_partial_windings().size() == 2) {
            sharedSectionFound = true;
        }
    }
    CHECK(sharedSectionFound);
    std::map<std::string, size_t> turnsPerWinding;
    for (auto& turn : coil.get_turns_description().value()) {
        turnsPerWinding[turn.get_winding()]++;
    }
    CHECK(turnsPerWinding["Primary"] == 12);
    CHECK(turnsPerWinding["Bias"] == 12);

    // Round trip: the serialized coil keeps the grouping.
    json serialized;
    to_json(serialized, coil);
    REQUIRE(serialized["functionalDescription"][0].contains("woundWith"));
    CHECK(serialized["functionalDescription"][0]["woundWith"] == json::array({"Bias"}));

    // The proportions/pattern overload — the one every web wind call uses —
    // must handle the grouped pair too (it virtualizes and devirtualizes).
    // The WASM wind path (wind_impl): default-constructed coil + setters,
    // wind-even-if-not-fit + include-additional-coordinates settings, then an
    // explicit compaction pass. Must survive a devirtualized grouped coil.
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_delimit_and_compact(true);
    settings.set_coil_include_additional_coordinates(true);
    OpenMagnetics::Coil coilLikeWasm;
    coilLikeWasm.set_bobbin_from_json(coilJson["bobbin"]);
    coilLikeWasm.set_functional_description(std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]));
    coilLikeWasm.preload_margins({});
    coilLikeWasm.set_layers_orientation(WindingOrientation::OVERLAPPING);
    coilLikeWasm.set_turns_alignment(CoilAlignment::CENTERED);
    REQUIRE(coilLikeWasm.wind(std::vector<double>{0.5, 0.5}, std::vector<size_t>{0, 1}, 1));
    REQUIRE(coilLikeWasm.delimit_and_compact());
    std::map<std::string, size_t> turnsPerWindingPattern;
    for (auto& turn : coilLikeWasm.get_turns_description().value()) {
        turnsPerWindingPattern[turn.get_winding()]++;
    }
    CHECK(turnsPerWindingPattern["Primary"] == 12);
    CHECK(turnsPerWindingPattern["Bias"] == 12);
}

// Single-window coils must NOT gain additionalCoordinates from the placement
// pass: StrayCapacitance/Temperature/CoilMesher interpret them as REAL
// conductor positions, so emitting the center-leg mirror would silently
// change their results for every classic design. The winding studio draws
// that mirror as pure display geometry instead (a regression here would be a
// physics change, not a display fix).
TEST_CASE("SingleWindow_CenterLegTurns_NoAdditionalCoordinates", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(false);
    settings.set_coil_include_additional_coordinates(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 10}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
    OpenMagnetics::Coil coil(coilJson, false);
    REQUIRE(coil.wind());
    auto woundTurns = coil.get_turns_description().value();
    REQUIRE(woundTurns.size() == 10);
    for (auto& turn : woundTurns) {
        CHECK(turn.get_coordinates()[0] > 0);
        CHECK(!turn.get_additional_coordinates());
    }
}

// Grouping through the web fixture (multi-window bobbin): reproduces the WASM
// wind path on the exact coil the studio sends.
TEST_CASE("WoundWith_WebFixtureMultiWindowBobbin", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_delimit_and_compact(true);
    settings.set_coil_include_additional_coordinates(true);
    std::ifstream fixtureFile("/home/alf/OpenMagnetics/WebFrontend/tests/fixtures/multicolumn_e42_transformer.json");
    REQUIRE(fixtureFile.is_open());
    json fixture = json::parse(fixtureFile);
    json coilJson = fixture["magnetic"]["coil"];
    coilJson.erase("sectionsDescription");
    coilJson.erase("layersDescription");
    coilJson.erase("turnsDescription");
    coilJson["functionalDescription"][0]["numberTurns"] = 12;
    coilJson["functionalDescription"][1]["numberTurns"] = 12;
    coilJson["functionalDescription"][1]["isolationSide"] = coilJson["functionalDescription"][0]["isolationSide"];
    coilJson["functionalDescription"][1]["wire"] = coilJson["functionalDescription"][0]["wire"];
    coilJson["functionalDescription"][0]["woundWith"] = json::array({"Secondary"});
    coilJson["functionalDescription"][1]["woundWith"] = json::array({"Primary"});

    OpenMagnetics::Coil coil;
    coil.set_bobbin_from_json(coilJson["bobbin"]);
    coil.set_functional_description(std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]));
    coil.preload_margins({});
    coil.set_layers_orientation(WindingOrientation::OVERLAPPING);
    coil.set_turns_alignment(CoilAlignment::CENTERED);
    REQUIRE(coil.wind(std::vector<double>{0.5, 0.5}, std::vector<size_t>{0, 1}, 1));
    REQUIRE(coil.delimit_and_compact());
    std::map<std::string, size_t> turnsPerWinding;
    for (auto& turn : coil.get_turns_description().value()) {
        turnsPerWinding[turn.get_winding()]++;
    }
    CHECK(turnsPerWinding["Primary"] == 12);
    CHECK(turnsPerWinding["Secondary"] == 12);
}

// Three windings, two sharing the main window plus one lateral: the shared
// window splits between the main-window windings and every turn lands.
TEST_CASE("MultiColumnWinding_ThreeWindings_SharedMainPlusLateral", "[constructive-model][coil][multi-column][smoke-test]") {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    auto columns = core.get_processed_description()->get_columns();
    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 20}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
    coilJson["functionalDescription"].push_back(json{{"name", "Tertiary"}, {"numberTurns", 8}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
    coilJson["functionalDescription"].push_back(json{{"name", "Secondary"}, {"numberTurns", 10}, {"numberParallels", 1},
                                                     {"isolationSide", "secondary"}, {"wire", "Round 0.475 - Grade 1"}});
    OpenMagnetics::Coil coil(coilJson, false);
    coil.get_mutable_functional_description()[2].set_winding_window(2);
    coil.set_core_columns(columns);
    REQUIRE(coil.wind());

    std::map<std::string, size_t> turnsPerWinding;
    auto threeWindingTurns = coil.get_turns_description().value();
    for (auto& turn : threeWindingTurns) {
        turnsPerWinding[turn.get_winding()]++;
        if (turn.get_winding() == "Secondary") {
            CHECK(turn.get_coordinates()[0] < 0);
        }
        else {
            CHECK(turn.get_coordinates()[0] > 0);
        }
    }
    CHECK(turnsPerWinding["Primary"] == 20);
    CHECK(turnsPerWinding["Tertiary"] == 8);
    CHECK(turnsPerWinding["Secondary"] == 10);

    // The two main-window sections split the window without overlapping.
    double primaryX = section_x(coil, "Primary");
    double tertiaryX = section_x(coil, "Tertiary");
    CHECK(primaryX > 0);
    CHECK(tertiaryX > 0);
    double primaryHalfWidth = 0;
    double tertiaryHalfWidth = 0;
    auto threeWindingSections = coil.get_sections_description().value();
    for (auto& section : threeWindingSections) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        if (section.get_partial_windings()[0].get_winding() == "Primary") {
            primaryHalfWidth = section.get_dimensions()[0] / 2;
        }
        if (section.get_partial_windings()[0].get_winding() == "Tertiary") {
            tertiaryHalfWidth = section.get_dimensions()[0] / 2;
        }
    }
    CHECK(std::abs(primaryX - tertiaryX) > primaryHalfWidth + tertiaryHalfWidth - 1e-9);
}
