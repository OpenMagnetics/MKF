// ABT #1290 (Alf, 2026-09-20): under REAL WINDING, MKF must not draw a conductor corner tighter
// than the wire's own minimum bend radius, and the limit must come from the wire's standard --
// "no magic numbers, all decisions are physics, manufacturing logic and standards".
//
// THESE TESTS ASSERT THE RULE, NOT A NUMBER. Nothing here is pinned to a particular radius: each
// case asks the WIRE itself (WireBend, i.e. IEC 60317-0-1 Table 6 read through IEC 60851-3
// 5.1.1) for its own flexibility radius, measures the standoff the COIL actually placed its
// innermost turn at, and then builds the two formers that put the drawn corner just below and
// just above that limit. A layout change moves both together, so the test keeps asking the same
// question; a test pinned to "0,215 mm throws" would pass for the wrong reason the day the
// layout moves.
//
// The former used is a RECTANGULAR column (an E core): a round or oblong column is all corner,
// its radius is its own, and the rule has nothing to judge there.
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "physical_models/WireBend.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <limits>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

// A former whose corner radius is exactly what the caller asks for, so the drawn corner can be
// placed either side of the wire's limit without touching anything else in the layout: the
// corner radius feeds the bend only -- the column's faces, the window and hence every turn's
// radial position stay where they were.
void set_former_corner_radius(OpenMagnetics::Coil& coil, double cornerRadius) {
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    processed.set_column_corner_radius(cornerRadius);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
}

// The radial gap between the former's face and the centre of the turn laid closest to it: the
// standoff the coil itself chose, read back off the wound geometry rather than assumed.
double innermost_turn_standoff(OpenMagnetics::Coil& coil) {
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    const auto processed = bobbin.get_processed_description().value();
    const double columnWidth = processed.get_column_width().value();
    double closest = std::numeric_limits<double>::max();
    // get_turns_description() hands the optional back BY VALUE, so the vector has to be a named
    // local or the range-for binds into a temporary that is already gone (-Werror=dangling-reference).
    const auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        closest = std::min(closest, std::abs(turn.get_coordinates()[0]));
    }
    return closest - columnWidth;
}

OpenMagnetics::Coil build_test_coil() {
    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.4 - Grade 1")});
    return OpenMagneticsTesting::get_quick_coil({6}, {1}, "E 25/13/7", 1,
                                               WindingOrientation::OVERLAPPING,
                                               WindingOrientation::OVERLAPPING,
                                               CoilAlignment::CENTERED,
                                               CoilAlignment::CENTERED, wires);
}

}  // namespace

TEST_CASE("Real winding refuses a corner tighter than the wire's own bend limit, and accepts one "
          "just above it",
          "[constructive-model][coil][real-geometry][wire-bend][abt1290]") {
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    auto coil = build_test_coil();
    bendSettings.set_coil_use_real_winding_geometry(true);

    // The wire's own limit, from the wire's own standard. Nothing in this test knows the number.
    auto wire = coil.resolve_wire(0);
    const double flexibilityRadius =
        WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND);
    REQUIRE(flexibilityRadius > 0);

    // A former comfortably above the limit: this must wind, and it is also how the standoff the
    // coil chooses is measured.
    set_former_corner_radius(coil, flexibilityRadius);
    REQUIRE_NOTHROW(coil.wind());
    REQUIRE(coil.get_turns_description());
    const double standoff = innermost_turn_standoff(coil);
    REQUIRE(standoff > 0);

    // The two formers either side of the limit. The drawn corner is the former's corner plus the
    // turn's standoff, so these put it at 0,99 and 1,01 times the wire's flexibility radius.
    const double justUnder = 0.99 * flexibilityRadius - standoff;
    const double justOver = 1.01 * flexibilityRadius - standoff;
    REQUIRE(justUnder > 0);  // both formers must be real formers, not negative radii
    REQUIRE(justOver > justUnder);

    SECTION("a corner at 0,99 of the wire's flexibility radius is refused") {
        set_former_corner_radius(coil, justUnder);
        REQUIRE_THROWS_AS(coil.wind(), UnwindableBendException);
    }

    SECTION("the refusal names the winding, the drawn radius, the limit and the criterion") {
        set_former_corner_radius(coil, justUnder);
        using Catch::Matchers::ContainsSubstring;
        REQUIRE_THROWS_WITH(coil.wind(),
                            ContainsSubstring(coil.get_functional_description()[0].get_name()) &&
                                ContainsSubstring("bend radius") &&
                                ContainsSubstring("FLEXIBILITY") &&
                                ContainsSubstring("IEC 60317-0-1"));
    }

    SECTION("a corner at 1,01 of the same radius is not refused") {
        set_former_corner_radius(coil, justOver);
        REQUIRE_NOTHROW(coil.wind());
        REQUIRE(coil.get_turns_description());
    }

    bendSettings.reset();
}

TEST_CASE("The ideal layout is not judged against the wire's bend limit",
          "[constructive-model][coil][wire-bend][abt1290]") {
    // Alf's scope, explicit: the rule is enforced ONLY behind the real-winding flag. The same
    // sharp former that is refused above must wind exactly as it always has with the flag off --
    // the classic 2D layout treats the column as a mathematically sharp rectangle by design, and
    // every existing result depends on that.
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    REQUIRE_FALSE(bendSettings.get_coil_use_real_winding_geometry());
    auto coil = build_test_coil();
    set_former_corner_radius(coil, 0.0);  // as sharp as a former can be asked to be
    REQUIRE_NOTHROW(coil.wind());
    REQUIRE(coil.get_turns_description());
}
