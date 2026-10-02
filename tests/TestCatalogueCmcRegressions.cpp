// WE common-mode chokes from the requirements-sheet catalogue that the winder or the stray-capacitance
// model refused when El Choker ran the whole catalogue. Fixtures cmc_catalogue_<part>.json are the
// enriched magnetics El Choker hands to magnetic_autocomplete (core, coil with stated sections, no turns).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <source_location>
#include <fstream>
#include "support/Settings.h"
#include "support/Utils.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/StrayCapacitance.h"
#include "TestingUtils.h"
#include "json.hpp"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;

namespace {
json load_cmc_catalogue_fixture(const std::string& partNumber) {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_catalogue_" + partNumber + ".json");
    std::ifstream file(path);
    REQUIRE(file.good());
    return json::parse(file);
}

// Reach of a toroidal turn's insulation from the ring axis, against the winding window's radial
// height (the bore surface the first ring rests on).
struct BoreReach {
    double excess;      // insulation reach minus the window's radial height; > 0 is inside the core
    size_t turnIndex;
};
BoreReach worst_bore_reach(OpenMagnetics::Coil& coil) {
    double radialHeight = coil.resolve_bobbin().get_processed_description()->get_winding_windows()[0].get_radial_height().value();
    auto turns = coil.get_turns_description().value();
    BoreReach worst{-std::numeric_limits<double>::infinity(), 0};
    for (size_t index = 0; index < turns.size(); ++index) {
        auto windingIndex = coil.get_winding_index_by_name(turns[index].get_winding());
        double reach = std::hypot(turns[index].get_coordinates()[0], turns[index].get_coordinates()[1]) +
                       coil.resolve_wire(windingIndex).get_maximum_outer_width() / 2;
        if (reach - radialHeight > worst.excess) {
            worst = {reach - radialHeight, index};
        }
    }
    return worst;
}

const std::vector<std::string> boreContactParts{"7448013501", "7448052502", "7448062603", "744841210"};
}

// WE common-mode chokes from the requirements-sheet catalogue: the winder rounded the radial
// height of the ring that rests on the bore to the NEAREST nanometre, which put its copper
// 0.12-0.40 nm past the window it was wound in, and the stray-capacitance model then refused the
// coil it had just been handed ("inside the core at the bore, 0.000000 mm past"). The ring now
// goes onto the grid on the clear side; nothing is forgiven in the check.
TEST_CASE("Test_Toroid_Bore_Ring_Is_Not_Rounded_Into_The_Core", "[constructive-model][coil][toroidal][cmc-catalogue][stray-capacitance]") {
    for (const auto& partNumber : boreContactParts) {
        INFO(partNumber);
        settings.reset();
        OpenMagnetics::Magnetic magnetic(load_cmc_catalogue_fixture(partNumber));
        auto result = magnetic_autocomplete(magnetic);
        REQUIRE(result.get_coil().get_turns_description());
        auto coil = result.get_coil();
        auto worst = worst_bore_reach(coil);
        INFO("worst excess " << worst.excess << " m, turn " << worst.turnIndex);
        CHECK(worst.excess <= 0);
        CHECK_NOTHROW(StrayCapacitance().calculate_capacitance(result));
    }
}

// The other side of the same check: a turn that really is inside the core is still refused. One
// nanometre is the winder's own grid, so an overlap of one grid step must not pass as rounding.
TEST_CASE("Test_Toroid_Turn_One_Nanometre_Inside_The_Bore_Is_Refused", "[constructive-model][coil][toroidal][cmc-catalogue][stray-capacitance]") {
    settings.reset();
    OpenMagnetics::Magnetic magnetic(load_cmc_catalogue_fixture(boreContactParts[0]));
    auto result = magnetic_autocomplete(magnetic);
    auto& coil = result.get_mutable_coil();
    auto worst = worst_bore_reach(coil);
    auto turns = coil.get_turns_description().value();
    auto coordinates = turns[worst.turnIndex].get_coordinates();
    double radius = std::hypot(coordinates[0], coordinates[1]);
    double pushedRadius = radius - worst.excess + 1e-9;
    turns[worst.turnIndex].set_coordinates({coordinates[0] * pushedRadius / radius, coordinates[1] * pushedRadius / radius});
    coil.set_turns_description(turns);
    REQUIRE(worst_bore_reach(coil).excess > 0.9e-9);
    REQUIRE_THROWS_WITH(StrayCapacitance().calculate_capacitance(result),
                        ContainsSubstring("lies inside the core at the bore"));
}

// WE 7448229004: 7 + 7 turns of 1.062 mm wire in the 3.7 mm bore of an epoxy-coated T14/8/9,
// one section per winding, each with a margin on both sides. The winder holds a section's margins
// as angles at its innermost radius. With 1.23 mm margins that leaves room for 4 + 3 turns in two
// rings; from 1.24 mm (the requirements sheet says 1.25) it leaves room for 6 in one ring and 6
// in two, and a third ring is narrower than the margin. The sizing used to walk ring counts one
// way, past that point, until the margin no longer fitted the ring and every ring was handed a
// forced turn; the coil came back unwound, blamed on "7 rings, 7.434 mm deep in a window of
// 3.7 mm". Unwound is the right verdict for that margin model, but the caller must be told the
// margins are what does not fit (ABT #930 contract: no turns, get_last_fit_failure says why).
TEST_CASE("Test_Toroid_Section_Margins_That_Do_Not_Fit_Are_Named", "[constructive-model][coil][toroidal][cmc-catalogue]") {
    auto with_margin = [](double margin) {
        auto magneticJson = load_cmc_catalogue_fixture("7448229004");
        for (auto& section : magneticJson["coil"]["sectionsDescription"]) {
            if (section.contains("margin")) {
                section["margin"] = {margin, margin};
            }
        }
        return OpenMagnetics::Magnetic(magneticJson);
    };
    SECTION("1.23 mm margins wind both windings in two rings") {
        settings.reset();
        auto result = magnetic_autocomplete(with_margin(1.23e-3));
        REQUIRE(result.get_coil().get_turns_description());
        CHECK(result.get_coil().get_turns_description()->size() == 14);
        CHECK(result.get_coil().get_last_fit_failure().empty());
    }
    for (double margin : {1.24e-3, 1.25e-3}) {
        DYNAMIC_SECTION("margins of " << margin * 1e3 << " mm do not fit, and the reason says so") {
            settings.reset();
            OpenMagnetics::Magnetic result;
            REQUIRE_NOTHROW(result = magnetic_autocomplete(with_margin(margin)));
            CHECK_FALSE(result.get_coil().get_turns_description());
            const auto reason = result.get_coil().get_last_fit_failure();
            INFO(reason);
            CHECK_THAT(reason, ContainsSubstring("winding 'L1' does not fit its round winding window: 7 turns of a 1.062 mm wire") &&
                               ContainsSubstring("margins of") && ContainsSubstring("innermost radius") &&
                               ContainsSubstring("with 1 ring they leave") && ContainsSubstring("with 2 rings they leave"));
            CHECK_THAT(reason, !ContainsSubstring("7 rings"));
        }
    }
}
