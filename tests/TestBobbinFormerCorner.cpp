// The corner a wire is wound round on a former (Alf, 2026-10-25, replacing the ABT #1303
// moulding-guide default "inside = 0.15 x wall"): the bobbin bore follows the core's own column, so
//     inside corner  = the core column's corner (from the core shape: a corner radius, EFD's
//                      chamfer, or 0 for a column MAS describes as a sharp rectangle)
//     outside corner = inside corner + wall thickness (Bayer, "Part and Mold Design", Fig. 2-4)
// and the bobbin throws when the corner cannot be known. No factor, no default.
//
// THESE TESTS ASSERT THE RULE, NOT A NUMBER: each reads the inside corner off the real core shape
// and the wall off the bobbin it built, and checks the corner the bobbin resolved against them.
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <algorithm>
#include <cmath>
#include <numbers>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

OpenMagnetics::Core processed_core(const std::string& shapeName) {
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, OpenMagneticsTesting::get_ground_gap(0.0005));
    if (!core.get_processed_description()) {
        core.process_data();
    }
    return core;
}

double shape_dimension(const OpenMagnetics::Core& core, const std::string& label) {
    auto dimensions = core.resolve_shape().get_dimensions().value();
    REQUIRE(dimensions.count(label) == 1);
    return resolve_dimensional_values(dimensions.at(label));
}

}  // namespace

TEST_CASE("A former on a square-cornered column (E) wraps its wire round the wall alone",
          "[constructive-model][bobbin][former-corner]") {
    // The E family dimensions only A-F: MAS describes its central column as a sharp rectangle.
    auto core = processed_core("E 25/13/7");
    const auto& centralColumn = core.get_processed_description()->get_columns()[0];
    REQUIRE(centralColumn.get_shape() == ColumnShape::RECTANGULAR);
    REQUIRE_FALSE(centralColumn.get_corner_radius());
    CHECK(OpenMagnetics::Bobbin::get_core_column_corner_radius(centralColumn, core, true) == 0);

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    const auto processed = bobbin.get_processed_description().value();
    const double wall = processed.get_column_thickness();
    // A real wall, or the outside corner would be 0 either way and the test vacuous.
    REQUIRE(wall > 0);
    // The datum is written when the bobbin is built, not invented when it is read.
    REQUIRE(processed.get_column_corner_radius());
    CHECK_THAT(bobbin.get_column_corner_radius(), Catch::Matchers::WithinAbs(wall, 1e-15));
}

TEST_CASE("A former on a chamfered column (EFD) follows the chamfer: inside corner q, outside q + wall",
          "[constructive-model][bobbin][former-corner]") {
    auto core = processed_core("EFD 25/13/9");
    const auto& centralColumn = core.get_processed_description()->get_columns()[0];
    REQUIRE(centralColumn.get_shape() == ColumnShape::IRREGULAR);
    REQUIRE_FALSE(centralColumn.get_corner_radius());
    const double chamfer = shape_dimension(core, "q");
    REQUIRE(chamfer > 0);
    CHECK(OpenMagnetics::Bobbin::get_core_column_corner_radius(centralColumn, core, true) == chamfer);

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    const double wall = bobbin.get_processed_description()->get_column_thickness();
    REQUIRE(wall > 0);
    CHECK_THAT(bobbin.get_column_corner_radius(), Catch::Matchers::WithinAbs(chamfer + wall, 1e-15));
}

TEST_CASE("A former on a column that carries its own corner (EPC, F2 / 2) carries it out by the wall",
          "[constructive-model][bobbin][former-corner]") {
    auto core = processed_core("EPC 17");
    const auto& centralColumn = core.get_processed_description()->get_columns()[0];
    REQUIRE(centralColumn.get_corner_radius());
    CHECK_THAT(centralColumn.get_corner_radius().value(),
               Catch::Matchers::WithinAbs(shape_dimension(core, "F2") / 2, 1e-9));

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    const double wall = bobbin.get_processed_description()->get_column_thickness();
    REQUIRE(wall > 0);
    CHECK_THAT(bobbin.get_column_corner_radius(),
               Catch::Matchers::WithinAbs(centralColumn.get_corner_radius().value() + wall, 1e-15));
}

TEST_CASE("An EFD shape without its chamfer dimension throws, naming the shape and 'q'",
          "[constructive-model][bobbin][former-corner]") {
    auto core = processed_core("EFD 25/13/9");
    auto shape = core.resolve_shape();
    auto dimensions = shape.get_dimensions().value();
    REQUIRE(dimensions.erase("q") == 1);
    shape.set_dimensions(dimensions);
    auto functionalDescription = core.get_functional_description();
    functionalDescription.set_shape(shape);
    // A fresh Core, so no shape resolved from the database before the edit is cached in it.
    MAS::MagneticCore withoutChamfer;
    withoutChamfer.set_functional_description(functionalDescription);
    withoutChamfer.set_processed_description(core.get_processed_description());
    OpenMagnetics::Core coreWithoutChamfer(withoutChamfer);
    REQUIRE(coreWithoutChamfer.resolve_shape().get_dimensions()->count("q") == 0);
    const auto& centralColumn = core.get_processed_description()->get_columns()[0];
    CHECK_THROWS_WITH(OpenMagnetics::Bobbin::get_core_column_corner_radius(centralColumn, coreWithoutChamfer, true),
                      Catch::Matchers::ContainsSubstring("EFD 25/13/9") &&
                      Catch::Matchers::ContainsSubstring("'q'"));
}

TEST_CASE("An irregular column the shape gives no corner for throws instead of defaulting",
          "[constructive-model][bobbin][former-corner]") {
    // ETD's lateral legs have a curved inner face: no corner radius describes them.
    auto core = processed_core("ETD 34/17/11");
    const auto& lateralColumn = core.get_processed_description()->get_columns()[1];
    REQUIRE(lateralColumn.get_shape() == ColumnShape::IRREGULAR);
    REQUIRE_FALSE(lateralColumn.get_corner_radius());
    CHECK_THROWS_WITH(OpenMagnetics::Bobbin::get_core_column_corner_radius(lateralColumn, core, false),
                      Catch::Matchers::ContainsSubstring("ETD 34/17/11"));
}

TEST_CASE("A rectangular former that states no corner and never met its core takes the interpolated "
          "corner until it is resolved against its core",
          "[constructive-model][bobbin][former-corner]") {
    auto core = processed_core("E 25/13/7");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    auto processed = bobbin.get_processed_description().value();
    processed.set_column_corner_radius(std::nullopt);
    bobbin.set_processed_description(processed);
    // No stated corner, no core shape named, no core: interpolated from the catalogue formers
    // around a rectangular core column, carried out by this former's own wall.
    CHECK_THAT(bobbin.get_column_corner_radius(),
               Catch::Matchers::WithinAbs(OpenMagnetics::Bobbin::get_interpolated_column_inside_corner_radius(
                                              ColumnShape::RECTANGULAR, processed.get_column_width().value()) +
                                              processed.get_column_thickness(),
                                          1e-15));
    // Resolved against its core, it gets the corner its bore follows.
    bobbin.set_column_corner_radius_from_core(core);
    CHECK_THAT(bobbin.get_column_corner_radius(),
               Catch::Matchers::WithinAbs(processed.get_column_thickness(), 1e-15));
}

TEST_CASE("A catalogue former that states no corner is given the chamfered column's one by autocomplete",
          "[constructive-model][bobbin][former-corner]") {
    auto core = processed_core("EFD 25/13/9");
    auto bobbin = find_bobbin_by_name("Bobbin EFD 25/13/9");
    REQUIRE(bobbin.get_processed_description());
    REQUIRE_FALSE(bobbin.get_processed_description()->get_column_corner_radius());
    const auto columnShape = bobbin.get_processed_description()->get_column_shape();
    REQUIRE((columnShape == ColumnShape::RECTANGULAR || columnShape == ColumnShape::IRREGULAR));
    bobbin.set_column_corner_radius_from_core(core);
    CHECK_THAT(bobbin.get_column_corner_radius(),
               Catch::Matchers::WithinAbs(shape_dimension(core, "q") +
                                              bobbin.get_processed_description()->get_column_thickness(),
                                          1e-15));
}

TEST_CASE("A former that declares its own corner radius is taken at its word",
          "[constructive-model][bobbin][former-corner]") {
    auto core = processed_core("E 25/13/7");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    auto processed = bobbin.get_processed_description().value();
    // Whatever the drawing says -- here, deliberately NOT what the rule would produce.
    const double declared = 3 * processed.get_column_thickness();
    processed.set_column_corner_radius(declared);
    bobbin.set_processed_description(processed);
    bobbin.set_column_corner_radius_from_core(core);
    CHECK(bobbin.get_column_corner_radius() == declared);
}

TEST_CASE("Under real winding, a turn on a square-cornered E former is charged the perimeter of a "
          "racetrack bent round the wall + its standoff",
          "[constructive-model][coil][real-geometry][former-corner]") {
    // The corner is only worth guarding where it reaches COPPER: this checks that the length MKF
    // charges each turn is the racetrack perimeter
    //     L = 4 d + 4 w + 8 s + (2 pi - 8) R,   R = (core column corner + wall) + s
    // (w, d the former's half-extents, s the turn's standoff off the former face), with the core
    // column corner of an E core 0. Every quantity on the right is read off the core, the bobbin
    // and the turn's own coordinate, not copied from a run. The removed 0.15 x wall default moves
    // each turn by (2 pi - 8) x 0.15 x wall, nine orders of magnitude above the tolerance.
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);

    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.4 - Grade 1")});
    auto coil = OpenMagneticsTesting::get_quick_coil({6}, {1}, "E 25/13/7", 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, wires);
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    const auto processed = bobbin.get_processed_description().value();
    REQUIRE(processed.get_column_shape() == ColumnShape::RECTANGULAR);
    const double wall = processed.get_column_thickness();
    REQUIRE(wall > 0);
    const double halfWidth = processed.get_column_width().value();
    const double halfDepth = processed.get_column_depth();

    auto core = processed_core("E 25/13/7");
    const double coreColumnCorner = OpenMagnetics::Bobbin::get_core_column_corner_radius(
        core.get_processed_description()->get_columns()[0], core, true);
    const double outsideCorner = coreColumnCorner + wall;

    REQUIRE_NOTHROW(coil.wind());
    const auto turns = coil.get_turns_description().value();
    size_t judged = 0;
    for (const auto& turn : turns) {
        // Real winding stores each winding's first turn as the lead-crossing station with length
        // 0 BY CONSTRUCTION (the realWindingCrossingBump in Coil::wind): it is not a wrap round
        // the former, so there is no racetrack to check. Every other turn is one.
        if (turn.get_length() == 0) {
            continue;
        }
        const double standoff = std::abs(turn.get_coordinates()[0]) - halfWidth;
        REQUIRE(standoff > 0);
        const double bendRadius = outsideCorner + standoff;
        REQUIRE(coil.get_turn_bend_radius(turn.get_name()));
        CHECK_THAT(coil.get_turn_bend_radius(turn.get_name()).value(),
                   Catch::Matchers::WithinAbs(bendRadius, 1e-12));
        const double expectedLength = 4 * halfDepth + 4 * halfWidth + 8 * standoff +
                                      (2 * std::numbers::pi - 8) * bendRadius;
        CHECK_THAT(turn.get_length(), Catch::Matchers::WithinAbs(expectedLength, 1e-12));
        ++judged;
    }
    // All but one crossing station per winding were checked, so this cannot pass vacuously.
    CHECK(judged + coil.get_functional_description().size() >= turns.size());
    REQUIRE(judged > 0);
    settings.reset();
}

namespace {

// A former as a coil-only caller hands it over: processed description only, no corner stated,
// no core shape named.
OpenMagnetics::Bobbin coreless_former(const std::string& shapeName) {
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(processed_core(shapeName), false);
    auto processed = bobbin.get_processed_description().value();
    processed.set_column_corner_radius(std::nullopt);
    OpenMagnetics::Bobbin coreless;
    coreless.set_processed_description(processed);
    REQUIRE_FALSE(coreless.get_functional_description());
    return coreless;
}

}  // namespace

TEST_CASE("The corner interpolator's samples are the catalogue formers' corners derived from the core "
          "shapes they name",
          "[constructive-model][bobbin][former-corner][corner-interpolator]") {
    const auto& samples = OpenMagnetics::Bobbin::get_column_corner_samples();
    REQUIRE(samples.count(ColumnShape::RECTANGULAR) == 1);
    REQUIRE(samples.count(ColumnShape::IRREGULAR) == 1);
    for (const auto& [columnShape, table] : samples) {
        json spelled;
        to_json(spelled, columnShape);
        size_t total = 0;
        for (const auto& [family, count] : table.samplesPerFamily) {
            INFO(spelled.get<std::string>() << " core column, bobbin family " << family << ": " << count << " samples");
            CHECK(count > 0);
            total += count;
        }
        INFO(spelled.get<std::string>() << ": " << total << " samples at " << table.columnWidths.size() << " widths");
        CHECK(total >= table.columnWidths.size());
        REQUIRE(table.columnWidths.size() == table.insideCornerRadii.size());
        REQUIRE_FALSE(table.columnWidths.empty());
        for (size_t i = 0; i + 1 < table.columnWidths.size(); ++i) {
            CHECK(table.columnWidths[i] < table.columnWidths[i + 1]);
        }
        // Piecewise linear: the knots are returned exactly, and a midpoint is the mean of its two.
        for (size_t i = 0; i < table.columnWidths.size(); ++i) {
            CHECK(OpenMagnetics::Bobbin::get_interpolated_column_inside_corner_radius(columnShape, table.columnWidths[i]) ==
                  table.insideCornerRadii[i]);
            if (i + 1 < table.columnWidths.size()) {
                const double midpoint = (table.columnWidths[i] + table.columnWidths[i + 1]) / 2;
                CHECK_THAT(OpenMagnetics::Bobbin::get_interpolated_column_inside_corner_radius(columnShape, midpoint),
                           Catch::Matchers::WithinAbs((table.insideCornerRadii[i] + table.insideCornerRadii[i + 1]) / 2, 1e-15));
            }
        }
    }

    // One sample read back against its source: the catalogue EFD 25/13/9 former's point is the q of
    // the EFD 25/13/9 shape at that former's column width (the mean, should other EFD formers share
    // the width).
    auto efdFormer = find_bobbin_by_name("Bobbin EFD 25/13/9");
    const double efdWidth = efdFormer.get_processed_description()->get_column_width().value();
    const auto& irregular = samples.at(ColumnShape::IRREGULAR);
    auto knot = std::find_if(irregular.columnWidths.begin(), irregular.columnWidths.end(),
                             [&](double width) { return std::abs(width - efdWidth) <= 1e-9; });
    REQUIRE(knot != irregular.columnWidths.end());
    REQUIRE(irregular.samplesPerFamily.count("efd") == 1);
    const double q = shape_dimension(processed_core("EFD 25/13/9"), "q");
    // Every EFD former of this width names an EFD shape; the point is q exactly unless another one
    // shares the width, which the knot count tells.
    if (irregular.samplesPerFamily.at("efd") == irregular.columnWidths.size()) {
        CHECK_THAT(irregular.insideCornerRadii[size_t(knot - irregular.columnWidths.begin())],
                   Catch::Matchers::WithinAbs(q, 1e-12));
    }
}

TEST_CASE("A former with no corner, no core shape and no core takes its corner interpolated from the "
          "catalogue formers around the same core column",
          "[constructive-model][bobbin][former-corner][corner-interpolator]") {
    // EFD's quick former has the core's IRREGULAR column: it interpolates among the EFD formers.
    auto bobbin = coreless_former("EFD 25/13/9");
    const auto processed = bobbin.get_processed_description().value();
    REQUIRE(processed.get_column_shape() == ColumnShape::IRREGULAR);
    const double interpolated = OpenMagnetics::Bobbin::get_interpolated_column_inside_corner_radius(
        ColumnShape::IRREGULAR, processed.get_column_width().value());
    // An EFD column has a chamfer: the interpolated inside corner is not zero.
    CHECK(interpolated > 0);
    CHECK_THAT(bobbin.get_column_corner_radius(),
               Catch::Matchers::WithinAbs(interpolated + processed.get_column_thickness(), 1e-15));
}

TEST_CASE("A coreless former whose column width is outside every catalogue span throws, naming it",
          "[constructive-model][bobbin][former-corner][corner-interpolator]") {
    auto bobbin = coreless_former("E 25/13/7");
    auto processed = bobbin.get_processed_description().value();
    processed.set_column_width(1.0);
    bobbin.set_processed_description(processed);
    CHECK_THROWS_WITH(bobbin.get_column_corner_radius(),
                      Catch::Matchers::ContainsSubstring("1000.000000 mm is outside") &&
                      Catch::Matchers::ContainsSubstring("rectangular"));
}

TEST_CASE("Corner precedence: a stated corner, then the core, then interpolation",
          "[constructive-model][bobbin][former-corner][corner-interpolator]") {
    auto core = processed_core("EFD 25/13/9");
    const double wall = OpenMagnetics::Bobbin::create_quick_bobbin(core, false).get_processed_description()->get_column_thickness();
    const double fromCore = shape_dimension(core, "q") + wall;

    SECTION("the core beats interpolation") {
        // EPC: its column carries a corner (F2 / 2), while every catalogue former around a
        // rectangular core column is an E (corner 0), so the two rules give different numbers.
        // (A quick former of a catalogue EFD shape sits on a catalogue knot and interpolates to
        // exactly its own q, which could not tell them apart.)
        auto epcCore = processed_core("EPC 17");
        auto bobbin = coreless_former("EPC 17");
        const auto processed = bobbin.get_processed_description().value();
        const double epcFromCore = epcCore.get_processed_description()->get_columns()[0].get_corner_radius().value() +
                                   processed.get_column_thickness();
        const double interpolated = bobbin.get_column_corner_radius();
        // Distinct, or this section could not tell the two apart.
        REQUIRE(std::abs(interpolated - epcFromCore) > 1e-6);
        bobbin.set_column_corner_radius_from_core(epcCore);
        CHECK_THAT(bobbin.get_column_corner_radius(), Catch::Matchers::WithinAbs(epcFromCore, 1e-15));
    }
    SECTION("a stated corner beats interpolation") {
        auto bobbin = coreless_former("EFD 25/13/9");
        auto processed = bobbin.get_processed_description().value();
        const double stated = 5 * wall;
        processed.set_column_corner_radius(stated);
        bobbin.set_processed_description(processed);
        CHECK(bobbin.get_column_corner_radius() == stated);
    }
    SECTION("the core shape a catalogue former names beats interpolation") {
        auto catalogueFormer = find_bobbin_by_name("Bobbin EFD 25/13/9");
        REQUIRE_FALSE(catalogueFormer.get_processed_description()->get_column_corner_radius());
        const auto processed = catalogueFormer.get_processed_description().value();
        const double byShape = catalogueFormer.get_column_corner_radius();
        CHECK_THAT(byShape, Catch::Matchers::WithinAbs(shape_dimension(core, "q") + processed.get_column_thickness(), 1e-15));
        // The same former stripped of its name falls to interpolation, over its own (rectangular)
        // processed column: a different number.
        OpenMagnetics::Bobbin unnamed;
        unnamed.set_processed_description(processed);
        CHECK(std::abs(unnamed.get_column_corner_radius() - byShape) > 1e-6);
    }
}

TEST_CASE("A magnetic built from MAS winds its coil after giving it the core, so a coreless former "
          "follows the core column and never falls to interpolation",
          "[constructive-model][magnetic][former-corner][corner-interpolator]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);

    // EPC 17: see the precedence test, its core corner and the interpolated one differ.
    auto core = processed_core("EPC 17");
    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.4 - Grade 1")});
    auto coil = OpenMagneticsTesting::get_quick_coil({6}, {1}, "EPC 17", 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, wires);
    auto bobbin = coreless_former("EPC 17");
    const auto processed = bobbin.get_processed_description().value();
    const double fromCore = core.get_processed_description()->get_columns()[0].get_corner_radius().value() +
                            processed.get_column_thickness();
    const double interpolated = bobbin.get_column_corner_radius();
    REQUIRE(std::abs(interpolated - fromCore) > 1e-6);

    json coilJson;
    to_json(coilJson, coil);
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson.erase("sectionsDescription");
    coilJson.erase("layersDescription");
    coilJson.erase("turnsDescription");
    coilJson.erase("groupsDescription");
    json coreJson;
    to_json(coreJson, core);
    json magneticJson = {{"core", coreJson}, {"coil", coilJson}};

    OpenMagnetics::Magnetic magnetic{magneticJson.get<MAS::Magnetic>()};
    const auto& woundCoil = magnetic.get_coil();
    REQUIRE(woundCoil.get_turns_description());
    const double halfWidth = processed.get_column_width().value();
    size_t judged = 0;
    const auto woundTurns = woundCoil.get_turns_description().value();
    for (const auto& turn : woundTurns) {
        if (turn.get_length() == 0) {
            continue;  // the lead-crossing station, see the E racetrack test above
        }
        auto bendRadius = woundCoil.get_turn_bend_radius(turn.get_name());
        REQUIRE(bendRadius);
        const double standoff = std::abs(turn.get_coordinates()[0]) - halfWidth;
        REQUIRE(standoff > 0);
        CHECK_THAT(bendRadius.value() - standoff, Catch::Matchers::WithinAbs(fromCore, 1e-12));
        ++judged;
    }
    REQUIRE(judged > 0);
    settings.reset();
}
