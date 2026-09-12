#include <source_location>
#include "constructive_models/Bobbin.h"
#include "constructive_models/Core.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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

// ABT #763: Bobbin::process_data() used to be non-deterministic on a bobbin json that
// carries no functionalDescription. BobbinDataProcessor::factory dereferenced the
// disengaged optional (`get_functional_description()->get_family()`), reading the
// optional's UNINITIALISED storage as a BobbinFamily: undefined behaviour whose outcome is
// decided by whatever the stack happened to hold. Through PyOpenMagnetics 1.7.0 the same
// unmodified catalogue row therefore threw the misleading "Unknown bobbin family" on 14 of
// 15 fresh processes, SEGFAULTED on others (reproducible on demand by padding the
// environment, which shifts the stack), and on the surviving run returned a
// 0.0 x 0.0 m winding window — a plausible zero, which is worse than the throw.
//
// The stack-scribbling below is what makes this catch the bug IN-PROCESS: a tight loop
// would otherwise keep re-reading the same benign leftover bytes and the UB would hide.
[[gnu::noinline]] void scribble_stack() {
    volatile unsigned char pattern[16384];
    for (size_t index = 0; index < sizeof(pattern); ++index) {
        pattern[index] = static_cast<unsigned char>(0xA5u + (index % 251u));
    }
}

TEST_CASE("process_data on a bobbin with no functionalDescription throws, deterministically (ABT #763)",
          "[constructive-model][bobbin][abt763]") {
    // Every one of these arrives at Bobbin(json) as an empty bobbin: get_stack_optional
    // finds nothing, so both descriptions stay disengaged. The first is the exact mistake
    // that surfaced this — a caller passing json.dumps(row), i.e. a json *string*, where a
    // json object is expected.
    const std::vector<json> notBobbinObjects = {
        json("{\"name\": \"Bobbin EER 48L horizontal 20-pin (Norwe N0037-186)\"}"),
        json::object(),
        json::array(),
        json(42),
        json({{"name", "no functional description here"}}),
    };

    SECTION("the refusal is stable across many stack states, and never returns a zero window") {
        for (const auto& notABobbin : notBobbinObjects) {
            for (int attempt = 0; attempt < 200; ++attempt) {
                scribble_stack();
                OpenMagnetics::Bobbin bobbin(notABobbin);
                REQUIRE_FALSE(bobbin.get_functional_description());
                // Pre-fix this threw "Unknown bobbin family" (most attempts), crashed, or
                // silently returned a 0 x 0 window. The message is the assertion: it is what
                // distinguishes "I read garbage" from "I know what is wrong with your input".
                CHECK_THROWS_WITH(bobbin.process_data(),
                                  Catch::Matchers::ContainsSubstring("no functionalDescription"));
                CHECK_THROWS_AS(bobbin.process_data(), OpenMagnetics::InvalidInputException);
                // Nothing may have been written: no plausible zero left behind.
                CHECK_FALSE(bobbin.get_processed_description());
            }
        }
    }

    SECTION("BobbinDataProcessor::factory refuses it too, not only the process_data wrapper") {
        for (int attempt = 0; attempt < 200; ++attempt) {
            scribble_stack();
            OpenMagnetics::Bobbin bobbin(json::object());
            CHECK_THROWS_WITH(OpenMagnetics::BobbinDataProcessor::factory(bobbin),
                              Catch::Matchers::ContainsSubstring("no functionalDescription"));
        }
    }

    SECTION("a real catalogue row processes to the SAME non-zero window on every evaluation") {
        // The determinism half: 200 fresh evaluations of one unmodified MAS row, each with a
        // freshly poisoned stack underneath, must agree bit for bit and must not be zero.
        auto reference = OpenMagnetics::find_bobbin_by_name("Bobbin ETD 49");
        REQUIRE(reference.get_functional_description());
        MAS::Bobbin referenceBase = reference;
        json referenceJson;
        to_json(referenceJson, referenceBase);
        referenceJson.erase("processedDescription");

        std::optional<double> firstWidth;
        std::optional<double> firstHeight;
        std::optional<double> firstArea;
        for (int attempt = 0; attempt < 200; ++attempt) {
            scribble_stack();
            OpenMagnetics::Bobbin bobbin(referenceJson);
            REQUIRE(bobbin.get_processed_description());
            auto windingWindows = bobbin.get_processed_description()->get_winding_windows();
            REQUIRE(windingWindows.size() == 1);
            REQUIRE(windingWindows[0].get_width());
            REQUIRE(windingWindows[0].get_height());
            REQUIRE(windingWindows[0].get_area());
            double width = windingWindows[0].get_width().value();
            double height = windingWindows[0].get_height().value();
            double area = windingWindows[0].get_area().value();
            CHECK(width > 0);
            CHECK(height > 0);
            CHECK(area > 0);
            if (!firstWidth) {
                firstWidth = width;
                firstHeight = height;
                firstArea = area;
            }
            INFO("attempt " << attempt);
            CHECK(width == firstWidth.value());
            CHECK(height == firstHeight.value());
            CHECK(area == firstArea.value());
        }
    }

    SECTION("a family whose processor finds none of its dimensions is an error, not a zero window") {
        // ABT #634's silent-zero mechanism, now loud: the E dimension set {e,f,s1,s2,l2}
        // read by the ETD processor, which asks for {d1,d2,d3,h1,h2}. flatten_dimensions
        // returns 0 for every missing key without complaining, so this used to produce a
        // perfectly plausible 0 x 0 winding window.
        auto eBobbin = OpenMagnetics::find_bobbin_by_name("Bobbin E13/4");
        REQUIRE(eBobbin.get_functional_description());
        REQUIRE(eBobbin.get_functional_description()->get_family() == BobbinFamily::E);
        auto functionalDescription = eBobbin.get_functional_description().value();
        functionalDescription.set_family(BobbinFamily::ETD);
        eBobbin.set_functional_description(functionalDescription);
        CHECK_THROWS_WITH(eBobbin.process_data(),
                          Catch::Matchers::ContainsSubstring("zero area"));
    }
}

}  // namespace

// ============================================================================
// ABT #1171 / WP2 — Bobbin::expand_pinout
//
// The catalogue's 342 pinouts were pin COUNTS: nothing downstream could place a
// pin, so no bobbin ever drew one. These lock the placement, the numbering and —
// just as importantly — the refusals: a count without a pitch, a pitch without
// pin dimensions, and an odd row asked to straddle a centralPitch all throw
// instead of producing a plausible-looking footprint nobody can trust.
// ============================================================================

namespace {

// The Miles-Platts PQ 20/16 footprint as data/bobbins.ndjson states it, and as the
// PC-B2016-14 datasheet draws it: 14 pins, 6 on one row and 8 on the other, rows
// 10.16 mm (0.4 in) apart, 3.81 mm and 2.54 mm pitches, 5.08 mm across the middle pair.
MAS::Pinout miles_platts_pq2016_pinout() {
    json pinout = json::parse(R"({
        "numberPins": 14,
        "numberPinsPerRow": [6, 8],
        "rowDistance": 0.01016,
        "pitch": [0.00381, 0.00254],
        "centralPitch": 0.00508,
        "pinDescription": {"shape": "round", "type": "tht", "dimensions": [0.00071, 0.00071, 0.00508]}
    })");
    MAS::Pinout parsed;
    MAS::from_json(pinout, parsed);
    return parsed;
}

MAS::WindingWindowElement test_winding_window(double height, double width) {
    MAS::WindingWindowElement windingWindow;
    windingWindow.set_height(height);
    windingWindow.set_width(width);
    windingWindow.set_area(height * width);
    windingWindow.set_coordinates(std::vector<double>({width / 2, 0}));
    return windingWindow;
}

double pin_x(const std::vector<MAS::Pin>& pins, const std::string& name) {
    for (const auto& pin : pins) {
        if (pin.get_name() && pin.get_name().value() == name) {
            return pin.get_coordinates().value()[0];
        }
    }
    throw std::runtime_error("no pin named " + name);
}

}  // namespace

TEST_CASE("expand_pinout places the Miles-Platts PQ 20/16 footprint (ABT #1171)",
          "[constructive-model][bobbin][pins][abt1171]") {
    // The record itself must still say what the test was written against.
    auto catalogueBobbin = OpenMagnetics::find_bobbin_by_name("Bobbin PQ 20/16");
    REQUIRE(catalogueBobbin.get_functional_description());
    REQUIRE(catalogueBobbin.get_functional_description()->get_pinout());
    auto cataloguePinout = catalogueBobbin.get_functional_description()->get_pinout().value();
    CHECK(cataloguePinout.get_number_pins() == 14);
    REQUIRE(cataloguePinout.get_number_pins_per_row());
    CHECK(cataloguePinout.get_number_pins_per_row().value() == std::vector<int64_t>({6, 8}));
    REQUIRE(cataloguePinout.get_row_distance());
    CHECK_THAT(cataloguePinout.get_row_distance().value(),
               Catch::Matchers::WithinAbs(0.01016, 1e-9));
    REQUIRE(cataloguePinout.get_central_pitch());
    CHECK_THAT(cataloguePinout.get_central_pitch().value(),
               Catch::Matchers::WithinAbs(0.00508, 1e-9));
    REQUIRE(cataloguePinout.get_pitch());
    REQUIRE(std::holds_alternative<std::vector<double>>(cataloguePinout.get_pitch().value()));
    auto cataloguePitches = std::get<std::vector<double>>(cataloguePinout.get_pitch().value());
    REQUIRE(cataloguePitches.size() == 2);
    CHECK_THAT(cataloguePitches[0], Catch::Matchers::WithinAbs(0.00381, 1e-9));
    CHECK_THAT(cataloguePitches[1], Catch::Matchers::WithinAbs(0.00254, 1e-9));

    // PQ 20/16: window 0.0116 high, wall 0.0007, column depth 0.0106 — the numbers only
    // move the pins along the column axis, so the row/pitch assertions do not depend on them.
    const double wallThickness = 0.0007;
    const double columnDepth = 0.0106;
    auto windingWindow = test_winding_window(0.0116, 0.0035);
    auto pins = OpenMagnetics::Bobbin::expand_pinout(miles_platts_pq2016_pinout(),
                                                     MAS::OrientationEnum::VERTICAL,
                                                     windingWindow, wallThickness, columnDepth);

    REQUIRE(pins.size() == 14);
    for (size_t index = 0; index < pins.size(); ++index) {
        REQUIRE(pins[index].get_name());
        CHECK(pins[index].get_name().value() == std::to_string(index + 1));
        REQUIRE(pins[index].get_coordinates());
        CHECK(pins[index].get_coordinates()->size() == 3);
        CHECK(pins[index].get_shape() == MAS::PinShape::ROUND);
        CHECK(pins[index].get_type() == MAS::PinDescriptionType::THT);
    }

    // Row 0 (6 pins, pitch 3.81, middle pair 5.08 apart) at z = -rowDistance/2, numbered
    // along +X; row 1 (8 pins, pitch 2.54) at +rowDistance/2, numbered back along -X.
    const std::vector<double> expectedRow0({-0.01016, -0.00635, -0.00254, 0.00254, 0.00635, 0.01016});
    for (size_t index = 0; index < expectedRow0.size(); ++index) {
        CHECK_THAT(pin_x(pins, std::to_string(index + 1)),
                   Catch::Matchers::WithinAbs(expectedRow0[index], 1e-9));
        CHECK_THAT(pins[index].get_coordinates().value()[2],
                   Catch::Matchers::WithinAbs(-0.00508, 1e-9));
    }
    const std::vector<double> expectedRow1({0.01016, 0.00762, 0.00508, 0.00254,
                                            -0.00254, -0.00508, -0.00762, -0.01016});
    for (size_t index = 0; index < expectedRow1.size(); ++index) {
        CHECK_THAT(pin_x(pins, std::to_string(index + 7)),
                   Catch::Matchers::WithinAbs(expectedRow1[index], 1e-9));
        CHECK_THAT(pins[6 + index].get_coordinates().value()[2],
                   Catch::Matchers::WithinAbs(0.00508, 1e-9));
    }

    // Vertical: the pins hang below the bottom flange's OUTER face, centred half a pin lower.
    const double expectedY = -(0.0116 / 2 + wallThickness + 0.00508 / 2);
    for (const auto& pin : pins) {
        CHECK_THAT(pin.get_coordinates().value()[1], Catch::Matchers::WithinAbs(expectedY, 1e-9));
        CHECK_FALSE(pin.get_rotation());
    }
}

TEST_CASE("expand_pinout refuses a pinout that is only a pin count (ABT #1171)",
          "[constructive-model][bobbin][pins][abt1171]") {
    auto windingWindow = test_winding_window(0.0116, 0.0035);

    SECTION("a pin count with no pitch and no row distance") {
        // 342 of the 504 catalogue records looked exactly like this before WP2.
        json onlyCount = json::parse(R"({"numberPins": 6})");
        MAS::Pinout pinout;
        MAS::from_json(onlyCount, pinout);
        CHECK_THROWS_AS(OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                             windingWindow, 0.0007, 0.0106),
                        OpenMagnetics::InvalidInputException);
        CHECK_THROWS_WITH(OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                               windingWindow, 0.0007, 0.0106),
                          Catch::Matchers::ContainsSubstring("no pitch"));
    }

    SECTION("a pitch with no row distance is still no footprint") {
        json noRowDistance = json::parse(R"({"numberPins": 6, "pitch": 0.00254})");
        MAS::Pinout pinout;
        MAS::from_json(noRowDistance, pinout);
        CHECK_THROWS_WITH(OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                               windingWindow, 0.0007, 0.0106),
                          Catch::Matchers::ContainsSubstring("no rowDistance"));
    }

    SECTION("geometry without pin dimensions cannot become a solid") {
        json noPinDescription = json::parse(
            R"({"numberPins": 12, "numberRows": 2, "rowDistance": 0.02032, "pitch": 0.00508})");
        MAS::Pinout pinout;
        MAS::from_json(noPinDescription, pinout);
        CHECK_THROWS_WITH(OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                               windingWindow, 0.0007, 0.0106),
                          Catch::Matchers::ContainsSubstring("no pinDescription"));
    }

    SECTION("an odd row cannot straddle a centralPitch") {
        json oddRow = json::parse(R"({
            "numberPins": 10, "numberPinsPerRow": [5, 5], "rowDistance": 0.02032,
            "pitch": 0.00508, "centralPitch": 0.00254,
            "pinDescription": {"shape": "round", "type": "tht", "dimensions": [0.0007, 0.0007, 0.004]}
        })");
        MAS::Pinout pinout;
        MAS::from_json(oddRow, pinout);
        CHECK_THROWS_WITH(OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                               windingWindow, 0.0007, 0.0106),
                          Catch::Matchers::ContainsSubstring("two MIDDLE pins"));
    }

    SECTION("a pin count that does not divide over its rows") {
        json odd = json::parse(R"({
            "numberPins": 13, "numberRows": 2, "rowDistance": 0.02032, "pitch": 0.00508,
            "pinDescription": {"shape": "round", "type": "tht", "dimensions": [0.0007, 0.0007, 0.004]}
        })");
        MAS::Pinout pinout;
        MAS::from_json(odd, pinout);
        CHECK_THROWS_WITH(OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                               windingWindow, 0.0007, 0.0106),
                          Catch::Matchers::ContainsSubstring("numberPinsPerRow"));
    }
}

TEST_CASE("expand_pinout matches the hand computation for a scraped ETD 29 pinout (ABT #1171)",
          "[constructive-model][bobbin][pins][abt1171]") {
    // The geometry the WP2 enrichment writes onto "Bobbin ETD 29 vertical 12-pin
    // (Norwe 90641-186)" from the Shulin row TF-2902: 12 pins, two rows, 20 mm apart,
    // 5 mm pitch. The scrape publishes no pin diameter or length, so the dimensions here
    // are the test's own, stated explicitly rather than defaulted anywhere in the code.
    json scraped = json::parse(R"({
        "numberPins": 12, "numberRows": 2, "rowDistance": 0.02, "pitch": 0.005,
        "pinDescription": {"shape": "round", "type": "tht", "dimensions": [0.0008, 0.0008, 0.0045]}
    })");
    MAS::Pinout pinout;
    MAS::from_json(scraped, pinout);

    // ETD 29/16/10 bobbin: window 0.0195 high and 0.00475 wide, wall 0.0008, column depth 0.0113.
    const double wallThickness = 0.0008;
    const double columnDepth = 0.0113;
    auto windingWindow = test_winding_window(0.0195, 0.00475);

    SECTION("vertical: pins under the bottom flange, rows across the depth") {
        auto pins = OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::VERTICAL,
                                                         windingWindow, wallThickness, columnDepth);
        REQUIRE(pins.size() == 12);
        // Six pins to a row, 5 mm apart, centred: -0.0125 .. +0.0125 by hand.
        const std::vector<double> expectedRow0({-0.0125, -0.0075, -0.0025, 0.0025, 0.0075, 0.0125});
        for (size_t index = 0; index < 6; ++index) {
            CHECK_THAT(pins[index].get_coordinates().value()[0],
                       Catch::Matchers::WithinAbs(expectedRow0[index], 1e-12));
            CHECK_THAT(pins[index].get_coordinates().value()[2],
                       Catch::Matchers::WithinAbs(-0.01, 1e-12));
            // -(0.0195/2 + 0.0008 + 0.0045/2)
            CHECK_THAT(pins[index].get_coordinates().value()[1],
                       Catch::Matchers::WithinAbs(-0.012800, 1e-12));
        }
        // Pin 7 turns the corner: same end of the board as pin 6, on the other row.
        CHECK_THAT(pins[6].get_coordinates().value()[0], Catch::Matchers::WithinAbs(0.0125, 1e-12));
        CHECK_THAT(pins[6].get_coordinates().value()[2], Catch::Matchers::WithinAbs(0.01, 1e-12));
        CHECK_THAT(pins[11].get_coordinates().value()[0], Catch::Matchers::WithinAbs(-0.0125, 1e-12));
        CHECK_THAT(pins[11].get_coordinates().value()[2], Catch::Matchers::WithinAbs(0.01, 1e-12));
    }

    SECTION("horizontal: pins off the two end flanges, rows along the column axis") {
        auto pins = OpenMagnetics::Bobbin::expand_pinout(pinout, MAS::OrientationEnum::HORIZONTAL,
                                                         windingWindow, wallThickness, columnDepth);
        REQUIRE(pins.size() == 12);
        for (size_t index = 0; index < 6; ++index) {
            CHECK_THAT(pins[index].get_coordinates().value()[1],
                       Catch::Matchers::WithinAbs(-0.01, 1e-12));
        }
        for (size_t index = 6; index < 12; ++index) {
            CHECK_THAT(pins[index].get_coordinates().value()[1],
                       Catch::Matchers::WithinAbs(0.01, 1e-12));
        }
        // -(0.0113/2 + 0.00475 + 0.0045/2) = -(0.00565 + 0.00475 + 0.00225)
        for (const auto& pin : pins) {
            CHECK_THAT(pin.get_coordinates().value()[2],
                       Catch::Matchers::WithinAbs(-0.01265, 1e-12));
            REQUIRE(pin.get_rotation());
            CHECK(pin.get_rotation().value() == std::vector<double>({90, 0, 0}));
        }
    }
}

TEST_CASE("A processed bobbin carries its pins and can be asked for one by name (ABT #1171)",
          "[constructive-model][bobbin][pins][abt1171]") {
    // A record with the full footprint AND the mounting orientation: the four things
    // process_data needs before it will place a pin.
    json bobbinJson = json::parse(R"({
        "name": "Bobbin WP2 Fixture",
        "functionalDescription": {
            "type": "standard", "family": "etd", "shape": "ETD 29/16/10",
            "orientation": "vertical",
            "dimensions": {"d1": {"nominal": 0.0212}, "d2": {"nominal": 0.0113},
                           "d3": {"nominal": 0.0102}, "h1": {"nominal": 0.0212},
                           "h2": {"nominal": 0.0195}, "s": {"nominal": 0.0008}},
            "pinout": {
                "numberPins": 12, "numberRows": 2, "rowDistance": 0.02032, "pitch": 0.00508,
                "pinDescription": {"shape": "round", "type": "tht",
                                   "dimensions": [0.0008, 0.0008, 0.0045]}
            }
        }
    })");
    OpenMagnetics::Bobbin bobbin(bobbinJson);
    REQUIRE(bobbin.get_processed_description());
    REQUIRE(bobbin.get_processed_description()->get_pins());
    CHECK(bobbin.get_processed_description()->get_pins()->size() == 12);

    auto pin = bobbin.get_pin("7");
    REQUIRE(pin.get_coordinates());
    CHECK_THAT(pin.get_coordinates().value()[2], Catch::Matchers::WithinAbs(0.01016, 1e-12));
    CHECK_THROWS_WITH(bobbin.get_pin("13"), Catch::Matchers::ContainsSubstring("no pin named '13'"));

    SECTION("the same record without an orientation places nothing rather than guessing") {
        json noOrientation = bobbinJson;
        noOrientation["functionalDescription"].erase("orientation");
        OpenMagnetics::Bobbin unoriented(noOrientation);
        REQUIRE(unoriented.get_processed_description());
        CHECK_FALSE(unoriented.get_processed_description()->get_pins());
        CHECK_THROWS_WITH(unoriented.get_pin("1"), Catch::Matchers::ContainsSubstring("has no pins"));
    }
}
