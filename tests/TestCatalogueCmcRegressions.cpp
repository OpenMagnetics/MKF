// WE common-mode chokes from the requirements-sheet catalogue that the winder or the stray-capacitance
// model refused when El Choker ran the whole catalogue. Fixtures cmc_catalogue_<part>.json are the
// enriched magnetics El Choker hands to magnetic_autocomplete (core, coil with stated sections, no turns).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <source_location>
#include <fstream>
#include <map>
#include <set>
#include <cmath>
#include <catch2/catch_approx.hpp>
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

// WE 7448229004: 7 + 7 turns of 1.062 mm wire in the bore of an epoxy-coated T14/8/9, one section
// per winding. The windings are separated by one plastic bar through the window, 2.5 mm thick (the
// requirements sheet and the datasheet drawing), so each section carries a 1.25 mm margin on both
// sides: the half of the bar on its side. The bar is a slab about the radial line at the section
// boundary, and a turn of radius r_w centred at radius rho clears it at an angle
// phi >= asin((t/2 + r_w) / rho) from that line.
namespace {
OpenMagnetics::Magnetic cmc_7448229004(std::optional<double> margin, bool coated) {
    auto magneticJson = load_cmc_catalogue_fixture("7448229004");
    if (margin) {
        for (auto& section : magneticJson["coil"]["sectionsDescription"]) {
            if (section.contains("margin")) {
                section["margin"] = {margin.value(), margin.value()};
            }
        }
    }
    if (!coated) {
        magneticJson["core"]["functionalDescription"].erase("coating");
        // The bobbin is the bore surface; without the coating that is the bare 4.0 mm bore of the T14/8/9.
        auto& window = magneticJson["coil"]["bobbin"]["processedDescription"]["windingWindows"][0];
        double boreRadius = OpenMagnetics::resolve_dimensional_values(
            magneticJson["core"]["functionalDescription"]["shape"]["dimensions"]["B"].get<MAS::DimensionWithTolerance>()) / 2;
        window["radialHeight"] = boreRadius;
        window["coordinates"][0] = boreRadius;
        window["area"] = std::numbers::pi * boreRadius * boreRadius;
    }
    return OpenMagnetics::Magnetic(magneticJson);
}

struct RingCount {
    std::map<std::string, std::set<double>> ringsPerWinding;
    std::map<std::string, std::map<double, size_t>> turnsPerRing;
};
RingCount count_rings(OpenMagnetics::Coil& coil) {
    RingCount count;
    const auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        // Ring by its distance from the axis, on a micrometre grid.
        double radius = std::round(std::hypot(turn.get_coordinates()[0], turn.get_coordinates()[1]) * 1e6) / 1e6;
        count.ringsPerWinding[turn.get_winding()].insert(radius);
        count.turnsPerRing[turn.get_winding()][radius]++;
    }
    return count;
}

}

TEST_CASE("Test_Toroid_Flat_Bar_Clearance_Angle_Is_Analytic", "[constructive-model][coil][toroidal][cmc-catalogue][flat-bar-margin]") {
    settings.reset();
    // The fixture's own wire, as the winder resolves it.
    auto wound = magnetic_autocomplete(cmc_7448229004(std::nullopt, true));
    auto wire = wound.get_mutable_coil().resolve_wire(0);
    double wireRadius = wire.get_maximum_outer_width() / 2;
    double halfBar = 1.25e-3;
    double ringRadius = 3.7e-3 - wireRadius;   // the first ring in the coated T14/8/9 bore
    double expected = asin((halfBar + wireRadius) / ringRadius) * 180 / std::numbers::pi;
    CHECK(OpenMagnetics::Coil::toroidal_bar_clearance_angle(halfBar, ringRadius, wire) == Catch::Approx(expected).epsilon(1e-12));
    // 34.20 deg for the 1.062 mm wire: rho sin(phi) is exactly the bar face plus the wire radius.
    CHECK(ringRadius * sin(expected * std::numbers::pi / 180) == Catch::Approx(halfBar + wireRadius).epsilon(1e-12));
    // The bar's angle on the ring is what it adds over a turn's own half pitch, asin(r_w / rho).
    CHECK(OpenMagnetics::Coil::toroidal_bar_margin_angle(halfBar, ringRadius, wire) ==
          Catch::Approx(expected - asin(wireRadius / ringRadius) * 180 / std::numbers::pi).epsilon(1e-12));
    CHECK(OpenMagnetics::Coil::toroidal_bar_margin_angle(0, ringRadius, wire) == 0);
    // An outer ring loses less angle than an inner one.
    CHECK(OpenMagnetics::Coil::toroidal_bar_margin_angle(halfBar, ringRadius, wire) <
          OpenMagnetics::Coil::toroidal_bar_margin_angle(halfBar, ringRadius - 2 * wireRadius, wire));
    // No room: the bar face plus the wire reach the turn's own radius.
    CHECK_THROWS_WITH(OpenMagnetics::Coil::toroidal_bar_clearance_angle(halfBar, halfBar + wireRadius, wire), ContainsSubstring("no room"));
}

// On the fixture's coated bore (3.7 mm radius) the first ring's turn centres sit at 3.169 mm; the
// bar faces leave 180 - 2 x 34.20 = 111.6 deg for the centres, and 7 turns need 6 pitches of
// 19.29 deg = 115.75 deg. So 6 fit in that ring, and the winding takes a second ring, whose centres
// at 2.107 mm hold 3 more. The turns all wind, in two rings.
TEST_CASE("Test_Toroid_Flat_Bar_7448229004_Coated_Winds_In_Two_Rings", "[constructive-model][coil][toroidal][cmc-catalogue][flat-bar-margin]") {
    settings.reset();
    auto result = magnetic_autocomplete(cmc_7448229004(std::nullopt, true));
    INFO(result.get_coil().get_last_fit_failure());
    REQUIRE(result.get_coil().get_turns_description());
    auto coil = result.get_coil();
    CHECK(coil.get_turns_description()->size() == 14);
    auto rings = count_rings(coil);
    for (const auto& [winding, radii] : rings.ringsPerWinding) {
        INFO(winding);
        CHECK(radii.size() == 2);
        for (const auto& [radius, turns] : rings.turnsPerRing[winding]) {
            UNSCOPED_INFO("FLATBAR " << winding << " ring at " << radius * 1e3 << " mm: " << turns << " turns");
        }
    }
    OpenMagneticsTesting::check_turns_clear_toroidal_bars(coil);
    CHECK(coil.get_last_fit_failure().empty());
}

// Without the coating the bore is 4.0 mm: the centres sit at 3.469 mm, the bar faces leave
// 118.2 deg and 7 turns need 105.7 deg, so each winding lies in a single ring with 12.6 deg spare.
TEST_CASE("Test_Toroid_Flat_Bar_7448229004_Uncoated_Winds_In_One_Ring", "[constructive-model][coil][toroidal][cmc-catalogue][flat-bar-margin]") {
    settings.reset();
    auto result = magnetic_autocomplete(cmc_7448229004(std::nullopt, false));
    INFO(result.get_coil().get_last_fit_failure());
    REQUIRE(result.get_coil().get_turns_description());
    auto coil = result.get_coil();
    CHECK(coil.get_turns_description()->size() == 14);
    auto rings = count_rings(coil);
    for (const auto& [winding, radii] : rings.ringsPerWinding) {
        INFO(winding);
        CHECK(radii.size() == 1);
    }
    OpenMagneticsTesting::check_turns_clear_toroidal_bars(coil);
}

// Margins wider than the bars of 7448229004 used to be "held as angles at the innermost radius",
// which charged every ring the deepest ring's angle; 1.24 and 1.25 mm came back unwound. As flat
// bars they wind (two rings in the coated bore), and so does 1.23 mm. What does not fit still says
// why, ring by ring: with 1.25 mm bars 10 turns need a third ring, whose centres at 1.045 mm are
// closer to the axis than the bar face plus the wire radius (1.781 mm), so there is no room.
TEST_CASE("Test_Toroid_Section_Margins_That_Do_Not_Fit_Are_Named", "[constructive-model][coil][toroidal][cmc-catalogue][flat-bar-margin]") {
    for (double margin : {1.23e-3, 1.24e-3, 1.25e-3}) {
        DYNAMIC_SECTION("flat bars of " << margin * 1e3 << " mm per side wind every turn") {
            settings.reset();
            auto result = magnetic_autocomplete(cmc_7448229004(margin, true));
            INFO(result.get_coil().get_last_fit_failure());
            REQUIRE(result.get_coil().get_turns_description());
            CHECK(result.get_coil().get_turns_description()->size() == 14);
            CHECK(result.get_coil().get_last_fit_failure().empty());
            auto coil = result.get_coil();
            OpenMagneticsTesting::check_turns_clear_toroidal_bars(coil);
        }
    }
    SECTION("10 turns per winding between 1.25 mm bars do not fit, and the reason says so") {
        settings.reset();
        auto magneticJson = load_cmc_catalogue_fixture("7448229004");
        for (auto& winding : magneticJson["coil"]["functionalDescription"]) {
            winding["numberTurns"] = 10;
        }
        OpenMagnetics::Magnetic result;
        REQUIRE_NOTHROW(result = magnetic_autocomplete(OpenMagnetics::Magnetic(magneticJson)));
        CHECK_FALSE(result.get_coil().get_turns_description());
        const auto reason = result.get_coil().get_last_fit_failure();
        INFO(reason);
        CHECK_THAT(reason, ContainsSubstring("does not fit its round winding window: 10 turns of a 1.062 mm wire") &&
                           ContainsSubstring("flat bars of 1.250 mm and 1.250 mm half-thickness") &&
                           ContainsSubstring("ring 1 (turn centres at 3.169 mm) holds 6") &&
                           ContainsSubstring("ring 2 (turn centres at 2.107 mm) holds 3") &&
                           ContainsSubstring("ring 3 (turn centres at 1.045 mm) has no room"));
    }
}
