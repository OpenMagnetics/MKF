// ABT #1303 (Alf, 2026-09-21): the corner a wire is wound round on a moulded former is the
// wall's OUTSIDE (convex) corner, and a wall of uniform thickness has
//     outside radius = inside radius + wall thickness
// (Bayer, "Part and Mold Design", p. 21, Fig. 2-4). The old fallback put the inside-corner rule
// (0.5 x wall) on the outside corner.
//
// THESE TESTS ASSERT THE RULE, NOT A NUMBER: each asks the bobbin for the corner it resolved and
// checks it against the wall it was moulded with and the inside corner from its source -- the
// core's own column corner when the core declares one, the moulding rule otherwise.
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <cmath>
#include <numbers>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

// A rectangular-column core whose MAS shape dimensions no column corner (the E family carries
// only A-F), so its bobbin has to resolve the corner rather than read it.
OpenMagnetics::Core rectangular_column_core() {
    auto core = OpenMagneticsTesting::get_quick_core("E 25/13/7", OpenMagneticsTesting::get_ground_gap(0.0005));
    if (!core.get_processed_description()) {
        core.process_data();
    }
    return core;
}

}  // namespace

TEST_CASE("A former with no declared corner wraps its wire round inside corner + wall, the inside "
          "corner from the moulding rule when the core gives none",
          "[constructive-model][bobbin][former-corner][abt1303]") {
    auto core = rectangular_column_core();
    const auto coreProcessed = core.get_processed_description().value();
    const auto& centralColumn = coreProcessed.get_columns()[0];
    REQUIRE(centralColumn.get_shape() == ColumnShape::RECTANGULAR);
    // The precondition that makes this the moulding-rule branch: the core says nothing.
    REQUIRE_FALSE(centralColumn.get_corner_radius());

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    const auto processed = bobbin.get_processed_description().value();
    REQUIRE_FALSE(processed.get_column_corner_radius());
    const double wall = processed.get_column_thickness();
    // A real wall, or the rule has nothing to act on and the test would pass vacuously.
    REQUIRE(wall > 0);

    const double outside = bobbin.get_column_corner_radius();
    const double inside = OpenMagnetics::Bobbin::get_moulded_inside_corner_radius(wall);
    CHECK_THAT(outside - wall, Catch::Matchers::WithinAbs(inside, 1e-12));
    // And the inside corner is the moulding rule's R/t, not some other fraction of the wall.
    CHECK_THAT(inside / wall, Catch::Matchers::WithinRel(OpenMagnetics::Bobbin::mouldedInsideCornerRadiusToWallThickness, 1e-12));
}

TEST_CASE("A former built on a core that declares its column corner carries that corner outwards "
          "by the wall",
          "[constructive-model][bobbin][former-corner][abt1303]") {
    auto core = rectangular_column_core();
    auto coreProcessed = core.get_processed_description().value();
    auto columns = coreProcessed.get_columns();
    // Any corner the column can carry: the rule, not this value, is under test.
    const double coreCorner = 0.25 * std::min(columns[0].get_width(), columns[0].get_depth());
    columns[0].set_corner_radius(coreCorner);
    coreProcessed.set_columns(columns);
    core.set_processed_description(coreProcessed);

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    const double wall = bobbin.get_processed_description()->get_column_thickness();
    REQUIRE(wall > 0);
    // The bore's inside corner IS the column's corner (the bore fits the column with no
    // clearance), and the outside corner the wire wraps is one wall further out.
    CHECK_THAT(bobbin.get_column_corner_radius(), Catch::Matchers::WithinAbs(coreCorner + wall, 1e-12));
}

TEST_CASE("A former that declares its own corner radius is taken at its word",
          "[constructive-model][bobbin][former-corner][abt1303]") {
    auto core = rectangular_column_core();
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    auto processed = bobbin.get_processed_description().value();
    // Whatever the drawing says -- here, deliberately NOT what either rule would produce.
    const double declared = 3 * processed.get_column_thickness();
    processed.set_column_corner_radius(declared);
    bobbin.set_processed_description(processed);
    CHECK(bobbin.get_column_corner_radius() == declared);
}

TEST_CASE("Under real winding, a turn on a rectangular former with no declared corner is charged "
          "the perimeter of a racetrack bent round inside corner + wall + its standoff",
          "[constructive-model][coil][real-geometry][former-corner][abt1303]") {
    // The fallback is only worth guarding where it reaches COPPER: this checks that the length
    // MKF charges each turn is the racetrack perimeter
    //     L = 4 d + 4 w + 8 s + (2 pi - 8) R,   R = outside corner + s
    // (w, d the former's half-extents, s the turn's standoff off the former face), and that the
    // outside corner in it is the moulded-wall rule, inside + wall. Every quantity on the right is
    // read off the bobbin and the turn's own coordinate, not copied from a run. Replacing the
    // fallback with the old 0.5 x wall moves each turn by (2 pi - 8) x 0.65 x wall (~1 mm here),
    // nine orders of magnitude above the tolerance.
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
    // Preconditions that make this the fallback path, on a real wall.
    REQUIRE(processed.get_column_shape() == ColumnShape::RECTANGULAR);
    REQUIRE_FALSE(processed.get_column_corner_radius());
    const double wall = processed.get_column_thickness();
    REQUIRE(wall > 0);
    const double halfWidth = processed.get_column_width().value();
    const double halfDepth = processed.get_column_depth();

    const double outsideCorner = OpenMagnetics::Bobbin::get_moulded_outside_corner_radius(
        OpenMagnetics::Bobbin::get_moulded_inside_corner_radius(wall), wall);

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
        // The radius the wind recorded as charged is this one: the fallback feeds the bend.
        REQUIRE(coil.get_turn_bend_radius(turn.get_name()));
        CHECK_THAT(coil.get_turn_bend_radius(turn.get_name()).value(),
                   Catch::Matchers::WithinAbs(bendRadius, 1e-12));
        // And the length is that racetrack's perimeter: the fallback feeds the copper.
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
