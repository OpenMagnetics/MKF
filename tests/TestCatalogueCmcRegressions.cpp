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
