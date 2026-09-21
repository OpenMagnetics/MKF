// ABT #1290 follow-up: MKF must RETAIN the bend radius it charged a turn's length for.
//
// WHY THIS EXISTS. Coil::get_turn_length_in_frame has always taken a turnBendRadius, and no
// caller ever passed one -- so the value_or default was always taken and MKF never recorded
// WHICH radius it charged. Nothing downstream, and nothing inside MKF, could find out. That is
// how MKF came to charge 15_gan's turns for a 0,665 mm corner while MVB++ drew a 0,219 mm one,
// with neither side able to see the disagreement.
//
// WHAT THESE TESTS ASSERT. Not a radius -- the IDENTITY between three things that used to be
// three separate re-derivations: the radius charged, the radius retained, and the length paid
// for. The length is reconstructed from the retained radius, so the assertions fail if the
// retained number is anything other than the one the length was built from.
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

OpenMagnetics::Coil build_bend_test_coil() {
    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.4 - Grade 1")});
    // A RECTANGULAR column: the only frame that has a corner distinct from the turn itself, and
    // therefore the only one for which MKF chooses a bend radius at all.
    return OpenMagneticsTesting::get_quick_coil({6}, {1}, "E 25/13/7", 1,
                                               WindingOrientation::OVERLAPPING,
                                               WindingOrientation::OVERLAPPING,
                                               CoilAlignment::CENTERED,
                                               CoilAlignment::CENTERED, wires);
}

// The racetrack length get_turn_length_in_frame charges under real winding, written out
// independently here so a retained radius that is not the charged one cannot reconstruct it.
double racetrack_length(const WoundColumnFrame& frame, double turnX, double bendRadius) {
    const double radius = frame.axisX == 0 ? turnX : std::abs(turnX - frame.axisX);
    const double standoff = radius - frame.columnWidth;
    return 4 * frame.columnDepth + 4 * frame.columnWidth + 8 * standoff +
           (2 * std::numbers::pi - 8) * bendRadius;
}

}  // namespace

TEST_CASE("The bend radius a turn's length was charged for is retained and reachable",
          "[constructive-model][coil][real-geometry][turn-bend-radius][abt1290]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto coil = build_bend_test_coil();
    settings.set_coil_use_real_winding_geometry(true);
    REQUIRE(coil.wind());
    REQUIRE(coil.get_turns_description());

    const auto turns = coil.get_turns_description().value();
    REQUIRE_FALSE(turns.empty());

    size_t reconstructed = 0;
    size_t zeroLengthStationsWithARadius = 0;
    for (const auto& turn : turns) {
        REQUIRE(turn.get_section());
        auto frame = coil.get_wound_column_frame_for_section(turn.get_section().value());
        // Only a cornered column charges a bend radius; nothing else must claim one.
        if (frame.shape != ColumnShape::RECTANGULAR && frame.shape != ColumnShape::IRREGULAR) {
            REQUIRE_FALSE(coil.get_turn_bend_radius(turn.get_name()));
            continue;
        }

        auto retained = coil.get_turn_bend_radius(turn.get_name());
        INFO("turn '" << turn.get_name() << "' x=" << turn.get_coordinates()[0]
             << " length=" << turn.get_length());
        REQUIRE(retained);

        // (a) the retained radius IS what the single derivation yields for this turn.
        const double turnX = std::abs(turn.get_coordinates()[0]);
        REQUIRE_THAT(retained.value(),
                     Catch::Matchers::WithinRel(coil.get_turn_bend_radius_in_frame(frame, turnX), 1e-12));

        if (turn.get_length() == 0) {
            // THE OPENING CROSSING OF A (layer, parallel). ABT #685/#674 deliberately zeroes its
            // LENGTH -- the copper between two layers is charged on the connection markers, and
            // leaving a full turn on the layer's opening station would count it twice. Its corner
            // was still charged and is still recorded, which is precisely why the radius has to be
            // retained rather than recovered from the length: from a zeroed station the length
            // says nothing about the bend.
            ++zeroLengthStationsWithARadius;
            continue;
        }

        // (b) the load-bearing one: the LENGTH this turn actually carries reconstructs from the
        // retained radius. If the coil retained any other number, this fails.
        REQUIRE_THAT(turn.get_length(),
                     Catch::Matchers::WithinRel(racetrack_length(frame, turnX, retained.value()), 1e-9));
        ++reconstructed;
    }
    // A silent zero-coverage pass is the failure mode these guard against: the test must have
    // seen both kinds of station, or it is not measuring what it claims.
    REQUIRE(reconstructed > 0);
    REQUIRE(zeroLengthStationsWithARadius > 0);

    settings.reset();
}

TEST_CASE("A caller's own solved bend radius is the one charged and the one handed back",
          "[constructive-model][coil][real-geometry][turn-bend-radius][abt1290]") {
    // THE value_or BRANCH. This is the path that was dead: a caller that has solved the turn's
    // real bend (WireBend, for a turn that lifts off its former) hands it in, and the length must
    // be charged for THAT radius, not for the former-corner-plus-standoff default -- and the call
    // must say which one it used.
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto coil = build_bend_test_coil();
    settings.set_coil_use_real_winding_geometry(true);
    REQUIRE(coil.wind());
    REQUIRE(coil.get_turns_description());

    const auto turns = coil.get_turns_description().value();
    const auto& turn = turns.front();
    REQUIRE(turn.get_section());
    auto frame = coil.get_wound_column_frame_for_section(turn.get_section().value());
    REQUIRE((frame.shape == ColumnShape::RECTANGULAR || frame.shape == ColumnShape::IRREGULAR));
    const double turnX = std::abs(turn.get_coordinates()[0]);

    const double defaultBend = coil.get_turn_bend_radius_in_frame(frame, turnX);
    REQUIRE(defaultBend > 0);
    // A radius a solver could plausibly return for a turn standing off its former, and far enough
    // from the default that the two lengths cannot be confused.
    const double solvedBend = 3 * defaultBend;

    // The derivation honours the override...
    REQUIRE_THAT(coil.get_turn_bend_radius_in_frame(frame, turnX, solvedBend),
                 Catch::Matchers::WithinRel(solvedBend, 1e-12));

    // ...the length is charged for it, and differs from the default length...
    double charged = std::numeric_limits<double>::quiet_NaN();
    auto overriddenLength = coil.get_turn_length_in_frame(frame, turnX, solvedBend, &charged);
    REQUIRE(overriddenLength);
    REQUIRE_THAT(overriddenLength.value(),
                 Catch::Matchers::WithinRel(racetrack_length(frame, turnX, solvedBend), 1e-9));

    auto defaultLength = coil.get_turn_length_in_frame(frame, turnX);
    REQUIRE(defaultLength);
    REQUIRE(std::abs(overriddenLength.value() - defaultLength.value()) > 1e-9);

    // ...and the call reports the radius it charged, which is the override, not the default.
    REQUIRE_FALSE(std::isnan(charged));
    REQUIRE_THAT(charged, Catch::Matchers::WithinRel(solvedBend, 1e-12));

    settings.reset();
}

TEST_CASE("Ideal winding charges no corner radius, so none is retained",
          "[constructive-model][coil][turn-bend-radius][abt1290]") {
    // Scope, unchanged from ABT #1290: the classic 2D layout treats the column as a mathematically
    // sharp rectangle and chooses no corner radius at all. Retaining one would be recording a
    // decision MKF never made.
    auto& settings = Settings::GetInstance();
    settings.reset();
    REQUIRE_FALSE(settings.get_coil_use_real_winding_geometry());
    auto coil = build_bend_test_coil();
    REQUIRE(coil.wind());
    REQUIRE(coil.get_turns_description());

    const auto turns = coil.get_turns_description().value();
    REQUIRE_FALSE(turns.empty());
    for (const auto& turn : turns) {
        REQUIRE_FALSE(coil.get_turn_bend_radius(turn.get_name()));
    }
}
