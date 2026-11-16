// A ONE-TURN WINDING IS AN OMEGA, AND A STORED LAYOUT'S LITZ IS COMPRESSED TO FIT IT.
//
// Omega (Alf, 2026-09-07: "a U or Omega symbol, just one crossing and the terminals, side by side, on the
// same height"; 2026-11-04: "U, legs side by side", centred on the face, "Contact is fine" for the omega
// legs only). Under real winding each parallel of a one-turn winding enters and leaves the same turn at
// the same station. When both its legs stand on one half of the face they are centred across that
// crossing, their two half extents apart, the entrance on the side its turn runs to (wrapDirection), so
// the turn never passes over its own exit leg. The two legs touch, and the parallels of one omega winding
// stand on the same two lanes at their own heights: neither is a meeting.
//
// Litz compression (Alf, 2026-11-03/04: "calculate the space it has available and calculate the litz
// compression in height for that ... and use the normal wind method"; MKF decides, no switch). When a
// stored layout's sections cannot hold the round litz bundle, the bundle is compressed in height to the
// largest k x OD at which every section of the winding holds it, keeping its area, and then fills the
// cell its wound layers leave across the section.
//
// Every expected value is derived here from the fixture and the wires, independently of the wind and of
// terminal_lane_plan.
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <optional>
#include <source_location>
#include <string>

using json = nlohmann::json;
using OpenMagnetics::ConnectionKind;
using OpenMagnetics::ConnectionRoute;
using MAS::TurnCrossSectionalShape;

namespace {

const std::filesystem::path here = std::filesystem::path{std::source_location::current().file_name()}.parent_path();

// omega_litz_er35_three_windings: ER 35/20/11 in N87 on MKF's basic bobbin, a stored stack of sections
// Primary / Secondary / Primary / Tertiary / Primary filling the 26.4 mm window but for 0.1 mm:
//   Primary   17 turns x 3 parallels of Litz 180x0.071 Grade 1 Single Served, one parallel per section,
//             each section 7.553333 mm tall;
//   Secondary 1 turn x 2 parallels of Rectangular 3x0.80 Grade 1, its section 1.77 mm tall (two flat
//             turns of 0.885 mm);
//   Tertiary  1 turn x 2 parallels of Rectangular 2.80x0.80 Grade 1, its section 1.77 mm tall.
json fixture() {
    std::ifstream file(here / "testData" / "omega_litz_er35_three_windings.json");
    REQUIRE(file.good());
    return json::parse(file);
}

OpenMagnetics::Magnetic wound_with_real_winding(const json& magneticJson) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    clear_databases();
    settings.set_coil_use_real_winding_geometry(true);
    return OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(magneticJson));
}

ConnectionRoute terminal_of(OpenMagnetics::Coil& coil, const std::string& winding, size_t parallel, ConnectionKind kind) {
    const auto routes = coil.get_connection_layout().routes;
    std::optional<ConnectionRoute> found;
    for (const auto& route : routes) {
        if (route.winding == winding && static_cast<size_t>(route.parallel) == parallel && route.kind == kind) {
            REQUIRE_FALSE(found);   // one entrance and one exit per parallel
            found = route;
        }
    }
    REQUIRE(found);
    return found.value();
}

}   // namespace

TEST_CASE("A stored layout whose sections cannot hold the round litz bundle compresses it, with no switch",
          "[constructive-model][coil][real-winding][litz][omega]") {
    const auto magneticJson = fixture()["magnetic"];
    double storedPrimaryHeight = 0;
    for (const auto& section : magneticJson["coil"]["sectionsDescription"]) {
        if (section["name"] == "Primary section 0") {
            storedPrimaryHeight = section["dimensions"][1].get<double>();
        }
    }
    REQUIRE(storedPrimaryHeight > 0);

    auto magnetic = wound_with_real_winding(magneticJson);
    REQUIRE_FALSE(OpenMagnetics::Settings::GetInstance().get_coil_litz_deformation());   // no switch
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());

    // The round bundle: OD 1.4235 mm. Each Primary section is the window's 6.0 mm wide, 4 bundles across,
    // and holds one parallel's 17 turns plus one crossing station per layer, 21 cross-sections in 4
    // layers. Every layer outside the first loses the slot its entrance lead runs out through, so 4
    // layers of R rows hold R + 3 (R - 1) = 4R - 3. Round, R = floor(7.5533 / 1.4235) = 5 rows hold 17:
    // it does not fit. 6 rows hold 21, exactly: the bundle is compressed to 7.5533 / 6 = 1.2589 mm tall,
    // k = 0.884, above the bundle's floor of 0.72, and it is 1.332 mm wide there, still 4 across.
    const auto primaryWire = OpenMagnetics::find_wire_by_name("Litz 180x0.071 - Grade 1 - Single Served");
    const double od = OpenMagnetics::Wire::round_litz_outer_diameter(primaryWire);
    const double height = storedPrimaryHeight / 6;
    REQUIRE(std::floor(storedPrimaryHeight / od) == 5);   // the premise: round, the section holds 5 rows
    REQUIRE(4 * 5 - 3 < 17 + 4);
    REQUIRE(4 * 6 - 3 == 17 + 4);
    REQUIRE(height / od >= OpenMagnetics::Wire::minimum_litz_flattening(primaryWire));
    REQUIRE(coil.get_litz_flattening(0));
    CHECK_THAT(coil.get_litz_flattening(0).value(), Catch::Matchers::WithinRel(height / od, 1e-9));

    // Its turns fill the cell their 4 layers leave across the section, keeping the bundle's area.
    const auto sections = coil.get_sections_description().value();
    std::optional<double> primaryWidth;
    for (const auto& section : sections) {
        if (section.get_name() == "Primary section 0") {
            primaryWidth = section.get_dimensions()[0];
        }
        // A one-turn winding is charged no closing station, so its stored height is all it needs and it
        // takes nothing from the Primary sections beside it.
        if (section.get_name() == "Secondary section 0" || section.get_name() == "Tertiary section 0") {
            INFO(section.get_name());
            CHECK_THAT(section.get_dimensions()[1], Catch::Matchers::WithinAbs(0.00177, 1e-9));
        }
    }
    REQUIRE(primaryWidth);
    const double width = primaryWidth.value() / 4;
    const auto turns = coil.get_turns_description().value();
    size_t primaryTurns = 0;
    for (const auto& turn : turns) {
        if (turn.get_winding() != "Primary") {
            continue;
        }
        ++primaryTurns;
        INFO(turn.get_name());
        REQUIRE(turn.get_cross_sectional_shape());
        CHECK(turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::RECTANGULAR);
        REQUIRE(turn.get_dimensions());
        CHECK_THAT(turn.get_dimensions().value()[0], Catch::Matchers::WithinAbs(width, 1e-8));   // the re-wind moves the section's width by under a nanometre
        CHECK_THAT(turn.get_dimensions().value()[1], Catch::Matchers::WithinRel(height, 1e-9));
    }
    CHECK(primaryTurns == 3 * (17 + 4));   // every parallel's 17 turns and its 4 crossing stations
    const double radius = OpenMagnetics::Wire::compressed_litz_corner_radius(od, width, height);
    CHECK_THAT(width * height - (4 - std::numbers::pi) * radius * radius,
               Catch::Matchers::WithinRel(std::numbers::pi * od * od / 4, 1e-9));
    OpenMagnetics::Settings::GetInstance().reset();
}

TEST_CASE("A one-turn winding's two legs stand centred across its crossing, the entrance on the side its turn runs to",
          "[constructive-model][coil][real-winding][omega]") {
    auto magnetic = wound_with_real_winding(fixture()["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto wires = coil.get_wires();
    size_t legsChecked = 0;
    for (const std::string winding : {"Secondary", "Tertiary"}) {
        const size_t windingIndex = coil.get_winding_index_by_name(winding);
        REQUIRE(coil.get_functional_description()[windingIndex].get_number_turns() == 1);   // the premise
        REQUIRE(coil.get_functional_description()[windingIndex].get_number_parallels() == 2);
        // A flat lead runs out flat: along the face it is the strip's coated width (3.085 and 2.885 mm).
        const double extent = wires[windingIndex].get_maximum_outer_width();
        std::optional<double> entranceLane, exitLane;
        for (size_t parallel = 0; parallel < 2; ++parallel) {
            const auto entrance = terminal_of(coil, winding, parallel, ConnectionKind::TERMINAL_ENTRANCE);
            const auto exit = terminal_of(coil, winding, parallel, ConnectionKind::TERMINAL_EXIT);
            INFO(winding << " parallel " << parallel);
            // One turn, one station: the entrance arrives where the exit leaves.
            REQUIRE_FALSE(entrance.toTurn.empty());
            REQUIRE(entrance.toTurn == exit.fromTurn);
            REQUIRE(entrance.side == exit.side);
            REQUIRE(entrance.waypoints.back() == exit.waypoints.front());
            REQUIRE(entrance.exitSlot);
            REQUIRE(exit.exitSlot);
            REQUIRE_FALSE(entrance.sleeveOuterDiameter);
            REQUIRE_FALSE(exit.sleeveOuterDiameter);
            REQUIRE(entrance.wrapDirection);
            const double towardsTheTurn = entrance.wrapDirection.value() > 0 ? 1.0 : -1.0;
            INFO("entrance lane " << entrance.exitSlot.value() * 1e3 << " mm, exit lane " << exit.exitSlot.value() * 1e3
                                  << " mm, wrap " << entrance.wrapDirection.value() << ", extent " << extent * 1e3 << " mm");
            // Centred across the crossing, each half its extent off it: the two legs touch.
            CHECK_THAT(entrance.exitSlot.value(), Catch::Matchers::WithinAbs(towardsTheTurn * extent / 2, 1e-12));
            CHECK_THAT(exit.exitSlot.value(), Catch::Matchers::WithinAbs(-towardsTheTurn * extent / 2, 1e-12));
            // Its parallels stand on the same two lanes, at their own heights.
            if (entranceLane) {
                CHECK(entrance.exitSlot.value() == entranceLane.value());
                CHECK(exit.exitSlot.value() == exitLane.value());
            }
            entranceLane = entrance.exitSlot.value();
            exitLane = exit.exitSlot.value();
            // Neither the touch nor the shared lanes is a meeting: no leg is reported.
            for (const auto* route : {&entrance, &exit}) {
                for (const auto& finding : route->findings) {
                    UNSCOPED_INFO(finding.witness);
                }
                CHECK(route->findings.empty());
                ++legsChecked;
            }
        }
    }
    CHECK(legsChecked == 8);
    OpenMagnetics::Settings::GetInstance().reset();
}

// The same stored stack with 18 Primary turns: 18 turns and 4 stations are 22 cross-sections. Unblocked,
// 4 layers of 6 rows (24) would hold them, and that is the compression a model blind to the leads' slots
// takes, k = 0.884; but each outer layer loses its lead's slot, so 4 x 6 - 3 = 21 do not. 7 rows (k = 0.758)
// make the bundle 1.53 mm wide, 3 across, and 3 layers of 7 rows less 2 blocked slots hold 19 < 21; 8 rows
// (k = 0.663) are below the bundle's floor of 0.72. No compression holds it: the round wind stands, and
// says why.
TEST_CASE("Where the leads' blocked slots leave no compression above its floor, the round litz wind stands and says why",
          "[constructive-model][coil][real-winding][litz][omega]") {
    auto magneticJson = fixture()["magnetic"];
    magneticJson["coil"]["functionalDescription"][0]["numberTurns"] = 18;
    auto magnetic = wound_with_real_winding(magneticJson);
    auto& coil = magnetic.get_mutable_coil();
    CHECK_FALSE(coil.get_litz_flattening(0));
    const auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        if (turn.get_winding() == "Primary") {
            INFO(turn.get_name());
            REQUIRE(turn.get_cross_sectional_shape());
            CHECK(turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND);
        }
    }
    OpenMagnetics::Settings::GetInstance().reset();
}
