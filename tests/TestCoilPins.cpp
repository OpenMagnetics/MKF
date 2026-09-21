// ABT #1172 (WP3, MAS RFC 0013): winding ends assigned to bobbin pins, the lead routed to its pin,
// and the manufacturability rules that read the assignment (R2, R3, R4, R14).
//
// The former is the real catalogue record "Bobbin PQ 32/30": vertical, 12 THT pins in two rows of
// 6 (pitch 5.08 mm, central pitch 7.62 mm), placed by Bobbin::expand_pinout from the pin rail
// (ABT #1207). It is one of the eight PQ records whose pinout MKF can place; no horizontal record
// can be placed today (ETD 59 has no pin-rail datum), so the horizontal convention is covered by
// the route and rule unit tests only. Vertical formers run row 0 from -X and row 1 from +X (the
// diagonal), so pin positions below are compared as WALK indexes: 0 at the row's start end.

#include "advisers/Manufacturability.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WireBend.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

const std::string formerBobbinName = "Bobbin PQ 32/30";
const std::string formerShapeName = "PQ 32/30";

struct Winding_ {
    std::string name;
    int64_t turns;
    int64_t parallels;
    std::string side;
    std::string wire;
    json connections = nullptr;
};

OpenMagnetics::Coil make_coil(const std::vector<Winding_>& windingSpecs,
                                    std::optional<OpenMagnetics::Inputs> inputs = std::nullopt,
                                    std::vector<size_t> pattern = {},
                                    size_t repetitions = 1) {
    json coilJson;
    coilJson["bobbin"] = formerBobbinName;
    coilJson["functionalDescription"] = json::array();
    for (const auto& spec : windingSpecs) {
        json winding;
        winding["name"] = spec.name;
        winding["numberTurns"] = spec.turns;
        winding["numberParallels"] = spec.parallels;
        winding["isolationSide"] = spec.side;
        winding["wire"] = spec.wire;
        if (!spec.connections.is_null()) {
            winding["connections"] = spec.connections;
        }
        coilJson["functionalDescription"].push_back(winding);
    }
    OpenMagnetics::Coil coil(coilJson, false);
    if (inputs) {
        coil.set_inputs(inputs.value());
    }
    if (pattern.empty()) {
        REQUIRE(coil.wind());
    }
    else {
        REQUIRE(coil.wind(pattern, repetitions));
    }
    return coil;
}

OpenMagnetics::Core former_core() {
    return OpenMagneticsTesting::get_quick_core(formerShapeName, json::parse("[]"), 1, "N87");
}

OpenMagnetics::Inputs reinforced_offline_inputs(double maximumVoltageRms = 400) {
    DimensionWithTolerance altitude;
    altitude.set_maximum(2000);
    DimensionWithTolerance mainSupplyVoltage;
    mainSupplyVoltage.set_nominal(230);
    return OpenMagneticsTesting::get_quick_insulation_inputs(altitude, Cti::GROUP_I, IsolationClass::REINFORCED, mainSupplyVoltage,
                                                             OvervoltageCategory::II, PollutionDegree::PD2,
                                                             {InsulationStandards::IEC_623681}, maximumVoltageRms, 1.5 * maximumVoltageRms, 100000);
}

ConnectionElement connection_of(const OpenMagnetics::Coil& coil, const std::string& winding, End end,
                                std::optional<int64_t> parallel = std::nullopt, size_t occurrence = 0) {
    for (const auto& w : coil.get_functional_description()) {
        if (w.get_name() != winding) {
            continue;
        }
        REQUIRE(w.get_connections());
        size_t seen = 0;
        const auto connections = w.get_connections().value();
        for (const auto& connection : connections) {
            if (connection.get_end() && connection.get_end().value() == end && connection.get_parallel() == parallel) {
                if (seen == occurrence) {
                    return connection;
                }
                ++seen;
            }
        }
    }
    FAIL("no such connection on winding " << winding);
    return ConnectionElement();
}

PlacedPin placed_pin(OpenMagnetics::Coil& coil, const std::string& name) {
    auto bobbin = coil.resolve_bobbin();
    for (const auto& pin : OpenMagnetics::Coil::place_pins(bobbin.get_processed_description()->get_pins().value())) {
        if (pin.name == name) {
            return pin;
        }
    }
    FAIL("no pin " << name);
    return PlacedPin();
}

// Position along the row counted from the row's START end: -X on row 0, +X on row 1 (vertical).
size_t walk(OpenMagnetics::Coil& coil, const std::string& name) {
    auto pin = placed_pin(coil, name);
    auto bobbin = coil.resolve_bobbin();
    size_t rowSize = 0;
    for (const auto& other : OpenMagnetics::Coil::place_pins(bobbin.get_processed_description()->get_pins().value())) {
        if (other.row == pin.row) {
            ++rowSize;
        }
    }
    return (!pin.hangsAlongZ && pin.row % 2 == 1) ? rowSize - 1 - pin.indexAlongRow : pin.indexAlongRow;
}

size_t walk_of(OpenMagnetics::Coil& coil, const std::string& winding, End end, std::optional<int64_t> parallel = std::nullopt) {
    return walk(coil, connection_of(coil, winding, end, parallel).get_pin_name().value());
}

}  // namespace

TEST_CASE("The PQ 32/30 former places 12 vertical pins in two rows of 6, and overlapping pins are refused (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto bobbin = find_bobbin_by_name(formerBobbinName);
    REQUIRE(bobbin.get_processed_description());
    REQUIRE(bobbin.get_processed_description()->get_pins());
    auto placed = OpenMagnetics::Coil::place_pins(bobbin.get_processed_description()->get_pins().value());
    REQUIRE(placed.size() == 12);
    CHECK_FALSE(placed.front().hangsAlongZ);
    size_t rowOne = 0;
    for (const auto& pin : placed) {
        rowOne += pin.row;
    }
    CHECK(rowOne == 6);

    // Two pins of one row whose centres are closer than their radii are not a footprint.
    auto pins = bobbin.get_processed_description()->get_pins().value();
    auto stacked = pins[1];
    stacked.set_name("stacked");
    auto coordinates = stacked.get_coordinates().value();
    coordinates[0] += 0.0002;
    stacked.set_coordinates(coordinates);
    pins.push_back(stacked);
    CHECK_THROWS_WITH(OpenMagnetics::Coil::place_pins(pins), Catch::Matchers::ContainsSubstring("overlap"));
}

TEST_CASE("Flyback on a vertical former: rows are isolation groups, start and finish adjacent, creepage met (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto inputs = reinforced_offline_inputs();
    auto coil = make_coil({{"Primary", 40, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Secondary", 6, 1, "secondary", "Round 0.5 - Grade 1"}},
                                inputs);
    auto core = former_core();
    auto result = coil.assign_pins(coil.resolve_bobbin(), core);
    CHECK_FALSE(result.skipped);

    auto primaryStart = placed_pin(coil, connection_of(coil, "Primary", End::START).get_pin_name().value());
    auto primaryFinish = placed_pin(coil, connection_of(coil, "Primary", End::FINISH).get_pin_name().value());
    auto secondaryStart = placed_pin(coil, connection_of(coil, "Secondary", End::START).get_pin_name().value());
    auto secondaryFinish = placed_pin(coil, connection_of(coil, "Secondary", End::FINISH).get_pin_name().value());

    // Rule 1: one row per side, the two sides on different rows.
    CHECK(primaryStart.row == primaryFinish.row);
    CHECK(secondaryStart.row == secondaryFinish.row);
    CHECK(primaryStart.row != secondaryStart.row);
    // Rule 2 and Wuerth's vertical convention: adjacent inner pins, starts in one diagonal corner.
    CHECK(walk_of(coil, "Primary", End::START) == 2);
    CHECK(walk_of(coil, "Primary", End::FINISH) == 3);
    CHECK(walk_of(coil, "Secondary", End::START) == 2);
    CHECK(walk_of(coil, "Secondary", End::FINISH) == 3);
    CHECK(primaryStart.centre[0] < primaryFinish.centre[0]);
    CHECK(secondaryStart.centre[0] > secondaryFinish.centre[0]);
    // Connection data from the pin.
    auto connection = connection_of(coil, "Primary", End::START);
    CHECK(connection.get_type().value() == ConnectionType::THT);
    CHECK_THAT(connection.get_diameter().value(), Catch::Matchers::WithinAbs(0.00091, 1e-9));
    CHECK_FALSE(connection.get_parallel());

    // Creepage: required from insulation coordination, met by the worst pin pair.
    REQUIRE(result.requiredCreepage);
    CHECK(result.requiredCreepage.value() >= 0.004);
    REQUIRE(result.worstPath);
    CHECK(result.worstPath->get_path() >= result.requiredCreepage.value());

    // Idempotent: a second run keeps every pin (they are now the design's own).
    auto before = coil.get_functional_description();
    coil.assign_pins(coil.resolve_bobbin(), core);
    CHECK(coil.has_complete_pin_connections());
    for (size_t index = 0; index < before.size(); ++index) {
        CHECK(before[index].get_connections()->size() == coil.get_functional_description()[index].get_connections()->size());
        for (size_t c = 0; c < before[index].get_connections()->size(); ++c) {
            CHECK(before[index].get_connections().value()[c].get_pin_name() ==
                  coil.get_functional_description()[index].get_connections().value()[c].get_pin_name());
        }
    }
}

TEST_CASE("A bobbin without pins is skipped and the connections stay untouched (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto coil = OpenMagneticsTesting::get_quick_coil({20, 5}, {1, 1}, "PQ 32/30");
    auto core = former_core();
    auto result = coil.assign_pins(coil.resolve_bobbin(), core);
    CHECK(result.skipped);
    CHECK_FALSE(coil.get_functional_description()[0].get_connections());
}

TEST_CASE("Creepage shortfall throws with the pin pair and the shortfall in mm (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    // Three sides on two rows: the third shares a row with the second, and a 1000 V rms reinforced
    // 62368-1 requirement (10 mm) cannot be met by adjacent pins 5.08 mm apart on one rail: the
    // secondary takes pins 7 and 8, the auxiliary 9 and 10, so 8 and 9 are neighbours.
    auto inputs = reinforced_offline_inputs(1000);
    auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1",
                                  json::parse(R"([{"end": "start", "pinName": "7"}, {"end": "finish", "pinName": "8"}])")},
                                 {"Auxiliary", 4, 1, "tertiary", "Round 0.5 - Grade 1",
                                  json::parse(R"([{"end": "start", "pinName": "9"}, {"end": "finish", "pinName": "10"}])")}},
                                inputs);
    auto core = former_core();
    CHECK_THROWS_WITH(coil.assign_pins(coil.resolve_bobbin(), core),
                      Catch::Matchers::ContainsSubstring("shortfall of") && Catch::Matchers::ContainsSubstring("mm") &&
                      Catch::Matchers::ContainsSubstring("'8'") && Catch::Matchers::ContainsSubstring("'9'"));
}

TEST_CASE("Three isolation sides on a two-row former split one row into blocks, deterministically (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto build = [&]() {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"},
                                     {"Auxiliary", 4, 1, "tertiary", "Round 0.5 - Grade 1"}});
        auto core = former_core();
        auto result = coil.assign_pins(coil.resolve_bobbin(), core);
        CHECK(result.creepageNotCheckedReason.find("inputs") != std::string::npos);
        return coil;
    };
    auto coil = build();
    auto primaryStart = placed_pin(coil, connection_of(coil, "Primary", End::START).get_pin_name().value());
    auto primaryFinish = placed_pin(coil, connection_of(coil, "Primary", End::FINISH).get_pin_name().value());
    auto secondaryStart = placed_pin(coil, connection_of(coil, "Secondary", End::START).get_pin_name().value());
    auto secondaryFinish = placed_pin(coil, connection_of(coil, "Secondary", End::FINISH).get_pin_name().value());
    auto auxiliaryStart = placed_pin(coil, connection_of(coil, "Auxiliary", End::START).get_pin_name().value());
    auto auxiliaryFinish = placed_pin(coil, connection_of(coil, "Auxiliary", End::FINISH).get_pin_name().value());

    // Sides in build order take the row with the fewest sides: primary row A, secondary row B,
    // tertiary back on row A. The secondary keeps its row and takes the inner pins; row A is
    // split into two blocks at the two rail ends, the primary (wound first) at the start end,
    // the two spare pins between them.
    CHECK(primaryStart.row != secondaryStart.row);
    CHECK(primaryStart.row == auxiliaryStart.row);
    CHECK(walk_of(coil, "Primary", End::START) == 0);
    CHECK(walk_of(coil, "Primary", End::FINISH) == 1);
    CHECK(walk_of(coil, "Auxiliary", End::START) == 4);
    CHECK(walk_of(coil, "Auxiliary", End::FINISH) == 5);
    CHECK(walk_of(coil, "Secondary", End::START) == 2);
    CHECK(walk_of(coil, "Secondary", End::FINISH) == 3);
    (void) primaryFinish;
    (void) secondaryFinish;
    (void) auxiliaryFinish;

    auto again = build();
    for (size_t index = 0; index < 3; ++index) {
        auto first = coil.get_functional_description()[index].get_connections().value();
        auto second = again.get_functional_description()[index].get_connections().value();
        REQUIRE(first.size() == second.size());
        for (size_t c = 0; c < first.size(); ++c) {
            CHECK(first[c].get_pin_name() == second[c].get_pin_name());
        }
    }
}

TEST_CASE("User pin names are honoured, and refused when they break the row plan or do not exist (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto core = former_core();
    SECTION("a valid user pin is kept, the rest is filled around it") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "3"}, {"end": "finish"}])")},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), core);
        CHECK(connection_of(coil, "Primary", End::START).get_pin_name().value() == "3");
        auto finish = placed_pin(coil, connection_of(coil, "Primary", End::FINISH).get_pin_name().value());
        CHECK(finish.row == placed_pin(coil, "3").row);
        auto secondary = placed_pin(coil, connection_of(coil, "Secondary", End::START).get_pin_name().value());
        CHECK(secondary.row != finish.row);
    }
    SECTION("a secondary pin on the primary's row throws") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "3"}, {"end": "finish"}])")},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "5"}, {"end": "finish"}])")}});
        CHECK_THROWS_WITH(coil.assign_pins(coil.resolve_bobbin(), core),
                          Catch::Matchers::ContainsSubstring("Wrong row") && Catch::Matchers::ContainsSubstring("pin '5'"));
    }
    SECTION("one side spread over both rows throws") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "3"}, {"end": "finish", "pinName": "10"}])")},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"}});
        CHECK_THROWS_WITH(coil.assign_pins(coil.resolve_bobbin(), core), Catch::Matchers::ContainsSubstring("rows 0 and 1"));
    }
    SECTION("a pin the bobbin does not have throws") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "99"}, {"end": "finish"}])")},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"}});
        CHECK_THROWS_WITH(coil.assign_pins(coil.resolve_bobbin(), core), Catch::Matchers::ContainsSubstring("'99'"));
    }
}

TEST_CASE("Bifilar windings share a pin while the wraps fit, else take adjacent pins with parallel set (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto core = former_core();
    SECTION("two 0.3 mm strands share their pins") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                     {"Secondary", 4, 2, "secondary", "Round 0.3 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), core);
        auto connections = coil.get_functional_description()[1].get_connections().value();
        CHECK(connections.size() == 2);
        for (const auto& connection : connections) {
            CHECK_FALSE(connection.get_parallel());
        }
    }
    SECTION("three 0.5 mm strands exceed two wires per pin and spread over adjacent pins") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                     {"Secondary", 4, 3, "secondary", "Round 0.5 - Grade 1"}});
        auto result = coil.assign_pins(coil.resolve_bobbin(), core);
        auto connections = coil.get_functional_description()[1].get_connections().value();
        CHECK(connections.size() == 6);
        std::vector<size_t> startIndexes;
        std::vector<size_t> finishIndexes;
        for (int64_t parallel = 0; parallel < 3; ++parallel) {
            startIndexes.push_back(walk_of(coil, "Secondary", End::START, parallel));
            finishIndexes.push_back(walk_of(coil, "Secondary", End::FINISH, parallel));
        }
        CHECK(startIndexes == std::vector<size_t>({0, 1, 2}));
        CHECK(finishIndexes == std::vector<size_t>({3, 4, 5}));
        bool noted = false;
        for (const auto& note : result.notes) {
            noted = noted || note.find("shorted on the PCB") != std::string::npos;
        }
        CHECK(noted);
    }
}

TEST_CASE("Series sections of one winding share a tap pin between start and finish (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    // Primary - Secondary - Primary: the primary is two series sections with one junction.
    auto coil = make_coil({{"Primary", 40, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Secondary", 6, 1, "secondary", "Round 0.5 - Grade 1"}},
                                std::nullopt, {0, 1}, 2);
    auto core = former_core();
    coil.assign_pins(coil.resolve_bobbin(), core);
    auto start = placed_pin(coil, connection_of(coil, "Primary", End::START).get_pin_name().value());
    auto tap = placed_pin(coil, connection_of(coil, "Primary", End::TAP).get_pin_name().value());
    auto finish = placed_pin(coil, connection_of(coil, "Primary", End::FINISH).get_pin_name().value());
    CHECK(start.row == tap.row);
    CHECK(tap.row == finish.row);
    CHECK(walk_of(coil, "Primary", End::TAP) == walk_of(coil, "Primary", End::START) + 1);
    CHECK(walk_of(coil, "Primary", End::FINISH) == walk_of(coil, "Primary", End::TAP) + 1);
    CHECK(coil.get_functional_description()[0].get_connections()->size() == 3);
}

TEST_CASE("The terminal lead is routed to its pin: routedLength grows by the run past the window (ABT #1172)",
          "[constructive-model][coil][pins][abt1172][winding-losses]") {
    Settings::GetInstance().set_coil_use_real_winding_geometry(true);
    auto coil = make_coil({{"Primary", 40, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Secondary", 6, 1, "secondary", "Round 0.5 - Grade 1"}});
    auto core = former_core();

    std::vector<ConnectionRoute> routesBefore;
    auto spacesBefore = coil.get_connection_reserved_spaces(&routesBefore);
    auto lengthsBefore = WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(coil);

    coil.assign_pins(coil.resolve_bobbin(), core);
    std::vector<ConnectionRoute> routesAfter;
    auto spacesAfter = coil.get_connection_reserved_spaces(&routesAfter);
    auto lengthsAfter = WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(coil);

    // Blocking is unchanged: the same rectangles, in the same place.
    REQUIRE(spacesBefore.size() == spacesAfter.size());
    for (size_t index = 0; index < spacesBefore.size(); ++index) {
        CHECK(spacesBefore[index].coordinates == spacesAfter[index].coordinates);
        CHECK(spacesBefore[index].dimensions == spacesAfter[index].dimensions);
    }
    REQUIRE(routesBefore.size() == routesAfter.size());
    size_t terminals = 0;
    for (size_t index = 0; index < routesAfter.size(); ++index) {
        const auto& before = routesBefore[index];
        const auto& after = routesAfter[index];
        bool terminal = after.kind == ConnectionKind::TERMINAL_ENTRANCE || after.kind == ConnectionKind::TERMINAL_EXIT;
        if (!terminal) {
            CHECK(after.pinName.empty());
            CHECK(after.routedLength == before.routedLength);
            continue;
        }
        ++terminals;
        INFO(after.winding << " " << (after.kind == ConnectionKind::TERMINAL_ENTRANCE ? "entrance" : "exit"));
        REQUIRE_FALSE(after.pinName.empty());
        CHECK(before.pinName.empty());
        CHECK(after.routedLength > before.routedLength);
        REQUIRE(after.pinWaypoints.size() >= 2);
        // The run ends on the axis of the assigned pin, at least a wire radius beyond the plane the
        // pin leaves (ABT #1237: the wire rests on the rail face), and starts at the window exit.
        auto pin = placed_pin(coil, after.pinName);
        const auto& pinEnd = after.kind == ConnectionKind::TERMINAL_ENTRANCE ? after.pinWaypoints.front() : after.pinWaypoints.back();
        CHECK_THAT(pinEnd[0], Catch::Matchers::WithinAbs(pin.base[0], 1e-12));
        CHECK(pinEnd[1] <= pin.base[1] - 0.00025);
        CHECK_THAT(pinEnd[2], Catch::Matchers::WithinAbs(pin.base[2], 1e-12));
        const auto& windowEnd = after.kind == ConnectionKind::TERMINAL_ENTRANCE ? after.pinWaypoints.back() : after.pinWaypoints.front();
        const auto& borderEnd = after.kind == ConnectionKind::TERMINAL_ENTRANCE ? after.waypoints.front() : after.waypoints.back();
        // The border point (radial, axial) leaves on the front face at its own axial level, at
        // least as far out as the border (plus the ride-over lift and its lane).
        CHECK_THAT(windowEnd[1], Catch::Matchers::WithinAbs(borderEnd[1], 1e-12));
        CHECK(windowEnd[2] <= -borderEnd[0] + 1e-12);
        double leg = 0;
        for (size_t k = 0; k + 1 < after.pinWaypoints.size(); ++k) {
            const auto& a = after.pinWaypoints[k];
            const auto& b = after.pinWaypoints[k + 1];
            leg += std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
        }
        CHECK_THAT(after.routedLength - before.routedLength, Catch::Matchers::WithinAbs(leg, 2e-9));
    }
    CHECK(terminals == 4);
    // The loss model reads the longer leads.
    for (size_t winding = 0; winding < 2; ++winding) {
        CHECK(lengthsAfter[winding][0] > lengthsBefore[winding][0]);
    }
}

TEST_CASE("route_leads_to_pins walks beside the pin, down to the rail face and into the pin (ABT #1172, #1237)",
          "[constructive-model][coil][pins][abt1172][abt1237]") {
    MAS::Pin vertical;
    vertical.set_name("1");
    vertical.set_shape(PinShape::ROUND);
    vertical.set_type(PinDescriptionType::THT);
    vertical.set_dimensions({0.001, 0.001, 0.004});
    vertical.set_coordinates(std::vector<double>({-0.005, -0.012, 0.006}));
    PinLeadRequest lead;
    lead.label = "lead";
    lead.pin = vertical;
    lead.windowExit = {0.010, 0.003};
    lead.diameter = 0.0005;
    // Sharp corners: this test is about the run's placement, not its bends (ABT #1296).
    lead.bendRadius = lead.diameter / 2;
    lead.exitX = 0.0;
    auto routes = OpenMagnetics::Coil::route_leads_to_pins({vertical}, {lead}, 0.0, 2);
    REQUIRE(routes.size() == 1);
    auto route = routes[0];
    // A one-pin row has no outside, so the lead comes from -X: beside the pin at
    // x = -0.005 - (0.0005 pin radius + 0.00025 wire radius + 0.0005 wire diameter) = -0.00625.
    // Exit (0, 0.003, -0.010) -> along x to -0.00625 (0.00625) -> down to one wire radius under the
    // pin base y = -0.010 (0.01325) -> across to the row z = 0.006 (0.016) -> into the pin (0.00125).
    REQUIRE(route.waypoints.size() == 5);
    CHECK_THAT(route.length, Catch::Matchers::WithinAbs(0.00625 + 0.01325 + 0.016 + 0.00125, 1e-12));
    CHECK_THAT(route.waypoints.front()[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(route.waypoints.front()[2], Catch::Matchers::WithinAbs(-0.010, 1e-12));
    CHECK_THAT(route.waypoints[1][0], Catch::Matchers::WithinAbs(-0.00625, 1e-12));
    CHECK_THAT(route.waypoints.back()[0], Catch::Matchers::WithinAbs(-0.005, 1e-12));
    CHECK_THAT(route.waypoints.back()[1], Catch::Matchers::WithinAbs(-0.01025, 1e-12));
    CHECK_THAT(route.waypoints.back()[2], Catch::Matchers::WithinAbs(0.006, 1e-12));

    // The ride-over lift moves the exit out along -Z.
    lead.lift = 0.001;
    route = OpenMagnetics::Coil::route_leads_to_pins({vertical}, {lead}, 0.0, 2)[0];
    CHECK_THAT(route.waypoints.front()[2], Catch::Matchers::WithinAbs(-0.011, 1e-12));
    lead.lift = 0;

    MAS::Pin horizontal = vertical;
    horizontal.set_coordinates(std::vector<double>({0.004, 0.010, -0.020}));
    horizontal.set_rotation(std::vector<double>({90, 0, 0}));
    lead.pin = horizontal;
    route = OpenMagnetics::Coil::route_leads_to_pins({horizontal}, {lead}, 0.0, 2)[0];
    // Along x to 0.004 - 0.00125 (0.00275), out along -Z to a wire radius beyond the base plane
    // z = -0.018 (0.00825), up to the row y = 0.010 (0.007), into the pin (0.00125).
    CHECK_THAT(route.length, Catch::Matchers::WithinAbs(0.00275 + 0.00825 + 0.007 + 0.00125, 1e-12));
    CHECK_THAT(route.waypoints.back()[2], Catch::Matchers::WithinAbs(-0.01825, 1e-12));
}

TEST_CASE("Autocomplete assigns pins when the bobbin has them (ABT #1172)", "[constructive-model][coil][pins][abt1172][masautocomplete]") {
    // ABT #1237: only with coil_connect_leads_to_pins on (off by default).
    SettingsGuard<bool> connectLeads(Settings::GetInstance(), &Settings::get_coil_connect_leads_to_pins,
                                     &Settings::set_coil_connect_leads_to_pins, true);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(former_core());
    json coilJson;
    coilJson["bobbin"] = formerBobbinName;
    coilJson["functionalDescription"] = json::parse(R"([
        {"name": "Primary", "numberTurns": 40, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 0.5 - Grade 1"},
        {"name": "Secondary", "numberTurns": 6, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 0.5 - Grade 1"}
    ])");
    magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
    auto completed = magnetic_autocomplete(magnetic);
    CHECK(completed.get_coil().has_complete_pin_connections());
    auto primary = completed.get_coil().get_functional_description()[0].get_connections();
    REQUIRE(primary);
    CHECK(primary->size() == 2);
}

// -----------------------------------------------------------------------------------------
// Manufacturability rules that read the assignment
// -----------------------------------------------------------------------------------------

namespace {
OpenMagnetics::Magnetic former_magnetic(OpenMagnetics::Coil coil) {
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(former_core());
    magnetic.set_coil(coil);
    return magnetic;
}
}  // namespace

TEST_CASE("R2: the default assignment passes, a finish at the start end and interleaved pins warn (ABT #1172)",
          "[adviser][manufacturability][pins][abt1172]") {
    Manufacturability manufacturability;
    SECTION("green: assign_pins' own plan") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                     {"Auxiliary", 6, 1, "primary", "Round 0.5 - Grade 1"},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), former_core());
        auto magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r2_pin_order(magnetic);
        INFO(finding.get_message());
        CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
    }
    SECTION("red: a start and finish swapped against the others") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "3"}, {"end": "finish", "pinName": "4"}])")},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "10"}, {"end": "finish", "pinName": "9"}])")}});
        auto magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r2_pin_order(magnetic);
        INFO(finding.get_message());
        // Row 0 is numbered along +X and row 1 back along -X, so on a vertical former the diagonal
        // convention reads 3 -> 4 and 9 -> 10 as the same direction; 10 -> 9 breaks it.
        CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("starts at the end of its rail"));
    }
    SECTION("red: interleaved pins of two windings on one rail") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "2"}, {"end": "finish", "pinName": "4"}])")},
                                     {"Auxiliary", 6, 1, "primary", "Round 0.5 - Grade 1",
                                      json::parse(R"([{"end": "start", "pinName": "3"}, {"end": "finish", "pinName": "5"}])")}});
        auto magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r2_pin_order(magnetic);
        INFO(finding.get_message());
        CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("cross"));
    }
    SECTION("not evaluated without pins") {
        auto coil = OpenMagneticsTesting::get_quick_coil({20, 5}, {1, 1}, "PQ 32/30");
        OpenMagnetics::Magnetic magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r2_pin_order(magnetic);
        CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
    }
}

TEST_CASE("R3: 28 and 34 AWG on one rail warn, 28 and 30 pass (ABT #1172)", "[adviser][manufacturability][pins][abt1172]") {
    Manufacturability manufacturability;
    auto rail = [&](const std::string& fineWire) {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.315 - Grade 1"},   // ~28 AWG
                                     {"Auxiliary", 6, 1, "primary", fineWire},
                                     {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), former_core());
        auto magnetic = former_magnetic(coil);
        return manufacturability.evaluate_r3_rail_gauge_spread(magnetic);
    };
    auto red = rail("Round 0.16 - Grade 1");    // ~34 AWG
    INFO(red.get_message());
    CHECK(red.get_status() == ManufacturabilityStatus::WARNING);
    CHECK(red.get_measured_value().value() > 5);
    auto green = rail("Round 0.25 - Grade 1");  // ~30 AWG
    INFO(green.get_message());
    CHECK(green.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("R4: wrap height against the standoff (ABT #1172)", "[adviser][manufacturability][pins][abt1172]") {
    Manufacturability manufacturability;
    auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Secondary", 4, 2, "secondary", "Round 0.5 - Grade 1"}});
    coil.assign_pins(coil.resolve_bobbin(), former_core());
    SECTION("not evaluated: the catalogue former states no standoff") {
        auto magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r4_wrap_height(magnetic);
        CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
        CHECK_THAT(finding.get_reason().value(), Catch::Matchers::ContainsSubstring("bobbin has no pin length/standoff"));
    }
    auto with_standoff = [&](double standoff) {
        auto bobbin = coil.resolve_bobbin();
        auto functionalDescription = bobbin.get_functional_description().value();
        BobbinBase base;
        base.set_mounting(OrientationEnum::HORIZONTAL);
        DimensionWithTolerance dimension;
        dimension.set_nominal(0.03);
        base.set_length(dimension);
        base.set_width(dimension);
        base.set_height(dimension);
        DimensionWithTolerance standoffDimension;
        standoffDimension.set_nominal(standoff);
        base.set_standoff(standoffDimension);
        functionalDescription.set_base(base);
        bobbin.set_functional_description(functionalDescription);
        auto withBase = coil;
        withBase.set_bobbin(bobbin);
        auto magnetic = former_magnetic(withBase);
        return manufacturability.evaluate_r4_wrap_height(magnetic);
    };
    // The secondary's two strands share a pin: 2 wires x 2 wrap turns x ~0.55 mm OD = ~2.2 mm.
    auto red = with_standoff(0.001);
    INFO(red.get_message());
    CHECK(red.get_status() == ManufacturabilityStatus::FAIL);
    CHECK(red.get_measured_value().value() > 0.002);
    auto green = with_standoff(0.003);
    INFO(green.get_message());
    CHECK(green.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("R14: shield ends on pins pass, a buried end warns, no shield is not applicable (ABT #1172)",
          "[adviser][manufacturability][pins][abt1172]") {
    Manufacturability manufacturability;
    auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Shield", 5, 1, "primary", "Round 0.5 - Grade 1"},
                                 {"Secondary", 4, 1, "secondary", "Round 0.5 - Grade 1"}});
    SECTION("no shield section: not applicable") {
        auto magnetic = former_magnetic(coil);
        CHECK(manufacturability.evaluate_r14_shield_terminations(magnetic).get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
    }
    // Declare the Shield winding's section as a shielding section.
    auto sections = coil.get_sections_description().value();
    for (auto& section : sections) {
        if (!section.get_partial_windings().empty() && section.get_partial_windings()[0].get_winding() == "Shield") {
            section.set_type(ElectricalType::SHIELDING);
        }
    }
    coil.set_sections_description(sections);
    SECTION("green: assign_pins puts both shield ends on pins") {
        auto result = coil.assign_pins(coil.resolve_bobbin(), former_core());
        auto magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r14_shield_terminations(magnetic);
        INFO(finding.get_message());
        CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
    }
    SECTION("red: the finish buried as a flying lead") {
        auto windings = coil.get_functional_description();
        for (auto& winding : windings) {
            if (winding.get_name() == "Shield") {
                winding.set_connections(json::parse(R"([{"end": "start", "pinName": "6"}, {"end": "finish", "type": "flyingLead"}])")
                                            .get<std::vector<ConnectionElement>>());
            }
        }
        coil.set_functional_description(windings);
        auto magnetic = former_magnetic(coil);
        auto finding = manufacturability.evaluate_r14_shield_terminations(magnetic);
        INFO(finding.get_message());
        CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("'Shield' finish"));
    }
}


// -----------------------------------------------------------------------------------------
// ABT #1237: the runs to the pins share no copper
// -----------------------------------------------------------------------------------------

namespace {

using Point3 = std::vector<double>;

double point_segment_distance(const Point3& p, const Point3& a, const Point3& b) {
    double ab[3], ap[3];
    double abab = 0, apab = 0;
    for (int k = 0; k < 3; ++k) {
        ab[k] = b[k] - a[k];
        ap[k] = p[k] - a[k];
        abab += ab[k] * ab[k];
        apab += ap[k] * ab[k];
    }
    double t = abab > 0 ? std::clamp(apab / abab, 0.0, 1.0) : 0.0;
    double d2 = 0;
    for (int k = 0; k < 3; ++k) {
        const double c = ap[k] - t * ab[k];
        d2 += c * c;
    }
    return std::sqrt(d2);
}

// Every point of `a`, sampled every 20 um along each segment, against every segment of `b`.
double sampled_distance(const std::vector<Point3>& a, const std::vector<Point3>& b, size_t skipLastOfA = 0) {
    double best = std::numeric_limits<double>::max();
    const size_t segments = a.size() < 2 + skipLastOfA ? 0 : a.size() - 1 - skipLastOfA;
    for (size_t i = 0; i < segments; ++i) {
        double length = 0;
        for (int k = 0; k < 3; ++k) {
            length += (a[i + 1][k] - a[i][k]) * (a[i + 1][k] - a[i][k]);
        }
        length = std::sqrt(length);
        const size_t samples = std::max<size_t>(2, size_t(std::ceil(length / 20e-6)) + 1);
        for (size_t s = 0; s < samples; ++s) {
            const double t = double(s) / double(samples - 1);
            Point3 p = {a[i][0] + t * (a[i + 1][0] - a[i][0]), a[i][1] + t * (a[i + 1][1] - a[i][1]), a[i][2] + t * (a[i + 1][2] - a[i][2])};
            for (size_t j = 0; j + 1 < b.size(); ++j) {
                best = std::min(best, point_segment_distance(p, b[j], b[j + 1]));
            }
        }
    }
    return best;
}

struct PinRun {
    std::string label;
    std::string pinName;
    std::vector<Point3> points;   // window exit -> pin axis
    double radius;
    double slotX;                 // Coil::terminal_exit_slots for this lead, recomputed from the routes
};

// ABT #1237: the exit slots by the documented rule, recomputed from the coil's routes, turns and wires.
std::vector<std::optional<double>> recompute_exit_slots(OpenMagnetics::Coil& coil, const std::vector<ConnectionRoute>& routes) {
    const auto wires = coil.get_wires();
    const auto turns = coil.get_turns_description().value();
    std::vector<double> diameters;
    std::vector<double> attachAxial;
    for (const auto& route : routes) {
        auto wire = wires[coil.get_winding_index_by_name(route.winding)];
        diameters.push_back(std::max({wire.get_maximum_outer_width(), wire.get_maximum_outer_height(), route.sleeveOuterDiameter.value_or(0.0)}));
        double axial = std::numeric_limits<double>::quiet_NaN();
        const std::string& turnName = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.toTurn : route.fromTurn;
        for (const auto& turn : turns) {
            if (turn.get_name() == turnName) {
                axial = turn.get_coordinates()[1];
            }
        }
        attachAxial.push_back(axial);
    }
    return OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attachAxial);
}

struct PinRunReport {
    double smallestRunSeparationOverRequired = std::numeric_limits<double>::max();
    double smallestRunSeparation = std::numeric_limits<double>::max();
    size_t runs = 0;
};

// Checks the ABT #1237 goals on a wound coil whose ends are assigned to pins, and returns the runs.
PinRunReport check_pin_runs_share_no_copper(OpenMagnetics::Coil& coil, std::vector<PinRun>* runsOut = nullptr) {
    std::vector<ConnectionRoute> routes;
    coil.get_connection_reserved_spaces(&routes);
    auto bobbin = coil.resolve_bobbin();
    const auto placedPins = OpenMagnetics::Coil::place_pins(bobbin.get_processed_description()->get_pins().value());
    const auto wires = coil.get_wires();
    const int64_t wrapTurns = OpenMagnetics::Coil::pin_wrap_turns();

    const auto slots = recompute_exit_slots(coil, routes);
    std::vector<PinRun> runs;
    for (size_t routeIndex = 0; routeIndex < routes.size(); ++routeIndex) {
        const auto& route = routes[routeIndex];
        if (route.pinName.empty()) {
            continue;
        }
        PinRun run;
        run.label = route.winding + " p" + std::to_string(route.parallel) + (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? " start" : " finish");
        run.pinName = route.pinName;
        run.points = route.pinWaypoints;
        if (route.kind == ConnectionKind::TERMINAL_ENTRANCE) {
            std::reverse(run.points.begin(), run.points.end());
        }
        auto wire = wires[coil.get_winding_index_by_name(route.winding)];
        run.radius = std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height()) / 2;
        REQUIRE(run.points.size() >= 3);
        REQUIRE(slots[routeIndex].has_value());
        run.slotX = slots[routeIndex].value();
        runs.push_back(run);
    }
    auto pin_named = [&](const std::string& name) -> const PlacedPin& {
        for (const auto& pin : placedPins) {
            if (pin.name == name) {
                return pin;
            }
        }
        FAIL("no pin " << name);
        return placedPins.front();
    };
    auto pin_axis = [](const PlacedPin& pin) {
        const size_t axis = pin.hangsAlongZ ? 2 : 1;
        auto tip = pin.centre;
        tip[axis] = 2 * pin.centre[axis] - pin.base[axis];
        return std::vector<Point3>{pin.base, tip};
    };
    // A wrap as the stretch of its pin's axis it covers, wrapTurns wire diameters from where the run arrives.
    auto wrap_axis = [&](const PinRun& run) {
        const auto& pin = pin_named(run.pinName);
        const size_t axis = pin.hangsAlongZ ? 2 : 1;
        auto top = run.points.back();
        auto bottom = top;
        bottom[axis] -= double(wrapTurns) * 2 * run.radius;
        return std::vector<Point3>{top, bottom};
    };

    PinRunReport report;
    report.runs = runs.size();
    for (size_t i = 0; i < runs.size(); ++i) {
        const auto& run = runs[i];
        INFO(run.label << " -> pin " << run.pinName);
        const auto& pin = pin_named(run.pinName);
        const size_t axis = pin.hangsAlongZ ? 2 : 1;
        // The run starts at MKF's exit slot, the x the in-window lead is drawn at, and leaves along -Z.
        CHECK_THAT(run.points[0][0], Catch::Matchers::WithinAbs(run.slotX, 1e-12));
        CHECK_THAT(run.points[1][0], Catch::Matchers::WithinAbs(run.slotX, 1e-12));
        // Ends on its own pin's axis, resting outside the rail face, on a level last leg.
        CHECK_THAT(run.points.back()[0], Catch::Matchers::WithinAbs(pin.centre[0], 1e-12));
        CHECK(run.points.back()[axis] <= pin.base[axis] - run.radius + 1e-12);
        CHECK_THAT(run.points[run.points.size() - 2][axis], Catch::Matchers::WithinAbs(run.points.back()[axis], 1e-12));
        // The wrap stays on the pin.
        CHECK(wrap_axis(run)[1][axis] >= pin_axis(pin)[1][axis] - 1e-12);
        // Never inside a pin rail (ABT #1249): sampled every 20 um.
        for (const auto& rail : bobbin.get_pin_rails()) {
            double closest = std::numeric_limits<double>::max();
            for (size_t k = 0; k + 1 < run.points.size(); ++k) {
                const double length = std::hypot(std::hypot(run.points[k + 1][0] - run.points[k][0], run.points[k + 1][1] - run.points[k][1]),
                                                 run.points[k + 1][2] - run.points[k][2]);
                const size_t samples = std::max<size_t>(2, size_t(std::ceil(length / 20e-6)) + 1);
                for (size_t sample = 0; sample < samples; ++sample) {
                    const double t = double(sample) / double(samples - 1);
                    double d2 = 0;
                    for (int c = 0; c < 3; ++c) {
                        const double p = run.points[k][c] + t * (run.points[k + 1][c] - run.points[k][c]);
                        const double excess = std::max(std::abs(p - rail.centre[c]) - rail.halfExtents[c], 0.0);
                        d2 += excess * excess;
                    }
                    closest = std::min(closest, std::sqrt(d2));
                }
            }
            INFO("rail " << rail.name);
            CHECK(closest >= run.radius - 1e-9);
        }
        // Clear of every foreign pin by its wrap radius plus a wire diameter; of its own pin up to the last leg.
        for (const auto& other : placedPins) {
            const bool own = other.name == run.pinName;
            const double found = sampled_distance(run.points, pin_axis(other), own ? 1 : 0);
            INFO("pin " << other.name);
            CHECK(found >= other.diameter / 2 + 3 * run.radius - 1e-9);
        }
        for (size_t j = i + 1; j < runs.size(); ++j) {
            const auto& otherRun = runs[j];
            INFO("against " << otherRun.label << " -> pin " << otherRun.pinName);
            const double required = run.radius + otherRun.radius;
            const double separation = std::min(sampled_distance(run.points, otherRun.points),
                                               sampled_distance(otherRun.points, run.points));
            CHECK(separation >= required - 1e-9);
            report.smallestRunSeparation = std::min(report.smallestRunSeparation, separation);
            report.smallestRunSeparationOverRequired = std::min(report.smallestRunSeparationOverRequired, separation / required);
            // Distinct exits, and drops (the leg along the pin axis) at distinct x.
            const double exitGap = std::hypot(std::hypot(run.points[0][0] - otherRun.points[0][0], run.points[0][1] - otherRun.points[0][1]),
                                              run.points[0][2] - otherRun.points[0][2]);
            CHECK(exitGap >= required - 1e-9);
            auto drop_x = [&](const PinRun& r) {
                for (size_t k = 0; k + 1 < r.points.size(); ++k) {
                    if (std::abs(r.points[k][0] - r.points[k + 1][0]) < 1e-12 &&
                        std::abs(r.points[k][axis == 1 ? 2 : 1] - r.points[k + 1][axis == 1 ? 2 : 1]) < 1e-12 &&
                        std::abs(r.points[k][axis] - r.points[k + 1][axis]) > 1e-12 && r.points[k + 1][axis] < r.points[k][axis]) {
                        return r.points[k][0];
                    }
                }
                FAIL("no drop along the pin axis in the run of " << r.label);
                return 0.0;
            };
            CHECK(std::abs(drop_x(run) - drop_x(otherRun)) >= required - 1e-9);
            // Runs against the other's wrap (a wrap is a band of radius pinRadius + wire radius).
            const auto& otherPin = pin_named(otherRun.pinName);
            const bool samePin = otherRun.pinName == run.pinName;
            CHECK(sampled_distance(run.points, wrap_axis(otherRun), samePin ? 1 : 0) - (otherPin.diameter / 2 + otherRun.radius) >= required - 1e-9);
            CHECK(sampled_distance(otherRun.points, wrap_axis(run), samePin ? 1 : 0) - (pin.diameter / 2 + run.radius) >= required - 1e-9);
        }
    }
    if (runsOut) {
        *runsOut = runs;
    }
    return report;
}

// ABT #1237: every run's exit x is its expected lane times its coated OD (one wire per side in these fixtures).
void check_exit_lanes(const std::vector<PinRun>& runs, const std::map<std::string, int>& expectedLanes) {
    CHECK(runs.size() == expectedLanes.size());
    for (const auto& run : runs) {
        INFO(run.label << " -> pin " << run.pinName);
        REQUIRE(expectedLanes.count(run.label));
        CHECK_THAT(run.points[0][0], Catch::Matchers::WithinAbs(double(expectedLanes.at(run.label)) * 2 * run.radius, 1e-12));
    }
}

json boost_inductor_magnetic_json() {
    // MVB++ tests/mas_complete_fixtures/boost_inductor_complete.json, the magnetic (ABT #1237's reporter).
    return json::parse(R"({
        "core": {"name": "PQ 20/16 N97 Gapped",
                 "functionalDescription": {"type": "twoPieceSet", "material": "N97", "shape": "PQ 26/25",
                                           "gapping": [{"type": "subtractive", "length": 0.0012},
                                                       {"type": "residual", "length": 1e-05},
                                                       {"type": "residual", "length": 1e-05}],
                                           "numberStacks": 1}},
        "coil": {"bobbin": "Bobbin PQ 26/25",
                 "functionalDescription": [{"name": "Primary", "numberTurns": 13, "numberParallels": 2,
                                            "isolationSide": "primary", "wire": "Litz 40x0.1 - Grade 1 - Single Served"}]}
    })");
}

}  // namespace

TEST_CASE("The boost inductor's four runs to PQ 26/25 pins share no copper (ABT #1237)",
          "[constructive-model][coil][pins][abt1237]") {
    SettingsGuard<bool> realWinding(Settings::GetInstance(), &Settings::get_coil_use_real_winding_geometry,
                                    &Settings::set_coil_use_real_winding_geometry, true);
    SettingsGuard<bool> connectLeads(Settings::GetInstance(), &Settings::get_coil_connect_leads_to_pins,
                                     &Settings::set_coil_connect_leads_to_pins, true);
    auto magneticJson = boost_inductor_magnetic_json();
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagnetics::Core(magneticJson.at("core")));
    magnetic.set_coil(OpenMagnetics::Coil(magneticJson.at("coil"), false));
    auto completed = magnetic_autocomplete(magnetic);
    auto coil = completed.get_coil();
    REQUIRE(coil.has_complete_pin_connections());
    // The assignment is unchanged by the routing: strands on adjacent inner pins, starts towards -X.
    CHECK(connection_of(coil, "Primary", End::START, 0).get_pin_name().value() == "2");
    CHECK(connection_of(coil, "Primary", End::START, 1).get_pin_name().value() == "3");
    CHECK(connection_of(coil, "Primary", End::FINISH, 0).get_pin_name().value() == "4");
    CHECK(connection_of(coil, "Primary", End::FINISH, 1).get_pin_name().value() == "5");

    std::vector<PinRun> runs;
    auto report = check_pin_runs_share_no_copper(coil, &runs);
    CHECK(report.runs == 4);
    CHECK(report.smallestRunSeparationOverRequired >= 1 - 1e-9);
    // The exit slots (Coil::terminal_exit_slots), coated OD 0.943 mm. The entrance bundle sits at the
    // plane side by side; the exit bundle cannot take the plane (each exit row is 4.7 um from its own
    // Z dragback's step-out row there), so it anchors one lane out. These are the attach slots MVB++'s
    // fan chose itself on MKF 10533d78 (entrances 0 / 0.943 mm, exits 0.943 / 1.886 mm).
    std::map<std::string, double> expectedSlot = {{"Primary p0 start", 0.0}, {"Primary p1 start", 0.943e-3},
                                                  {"Primary p0 finish", 0.943e-3}, {"Primary p1 finish", 1.886e-3}};
    for (const auto& run : runs) {
        INFO(run.label);
        REQUIRE(expectedSlot.count(run.label));
        // 0.943 mm is the served litz's coated OD to the micrometre; the lane arithmetic is exact.
        CHECK_THAT(run.points[0][0], Catch::Matchers::WithinAbs(expectedSlot.at(run.label), 1e-7));
        CHECK_THAT(run.points[0][0] / (2 * run.radius), Catch::Matchers::WithinAbs(std::round(expectedSlot.at(run.label) / 0.943e-3), 1e-9));
    }
}

TEST_CASE("The PQ 32/30 flyback's runs share no copper, on both rows and with strands sharing a pin (ABT #1237)",
          "[constructive-model][coil][pins][abt1237]") {
    SettingsGuard<bool> realWinding(Settings::GetInstance(), &Settings::get_coil_use_real_winding_geometry,
                                    &Settings::set_coil_use_real_winding_geometry, true);
    auto core = former_core();
    SECTION("primary and secondary on opposite rows") {
        auto coil = make_coil({{"Primary", 40, 1, "primary", "Round 0.5 - Grade 1"},
                               {"Secondary", 6, 1, "secondary", "Round 0.5 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), core);
        std::vector<PinRun> runs;
        auto report = check_pin_runs_share_no_copper(coil, &runs);
        CHECK(report.runs == 4);
        // Exit lanes (x / coated OD): the Primary's finish leaves its own dragback's plane one lane out.
        check_exit_lanes(runs, {{"Primary p0 start", 0}, {"Primary p0 finish", 1}, {"Secondary p0 start", 0}, {"Secondary p0 finish", 0}});
    }
    SECTION("two strands share each pin") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                               {"Secondary", 4, 2, "secondary", "Round 0.3 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), core);
        CHECK_FALSE(connection_of(coil, "Secondary", End::START).get_parallel().has_value());
        std::vector<PinRun> runs;
        auto report = check_pin_runs_share_no_copper(coil, &runs);
        CHECK(report.runs == 6);
        check_exit_lanes(runs, {{"Primary p0 start", 0}, {"Primary p0 finish", 0}, {"Secondary p0 start", 0}, {"Secondary p1 start", 1},
                                {"Secondary p0 finish", 0}, {"Secondary p1 finish", 1}});
    }
    SECTION("three strands on adjacent pins") {
        auto coil = make_coil({{"Primary", 20, 1, "primary", "Round 0.5 - Grade 1"},
                               {"Secondary", 4, 3, "secondary", "Round 0.5 - Grade 1"}});
        coil.assign_pins(coil.resolve_bobbin(), core);
        std::vector<PinRun> runs;
        auto report = check_pin_runs_share_no_copper(coil, &runs);
        CHECK(report.runs == 8);
        check_exit_lanes(runs, {{"Primary p0 start", 0}, {"Primary p0 finish", 0}, {"Secondary p0 start", 0}, {"Secondary p1 start", 1},
                                {"Secondary p2 start", 2}, {"Secondary p0 finish", 0}, {"Secondary p1 finish", 1}, {"Secondary p2 finish", 2}});
    }
}

TEST_CASE("A lead that cannot be routed clear throws naming what it collides with (ABT #1237)",
          "[constructive-model][coil][pins][abt1237]") {
    auto make_pin = [](const std::string& name, double x, double length) {
        MAS::Pin pin;
        pin.set_name(name);
        pin.set_shape(PinShape::ROUND);
        pin.set_type(PinDescriptionType::THT);
        pin.set_dimensions({0.0008, 0.0008, length});
        pin.set_coordinates(std::vector<double>({x, -0.012 - length / 2, -0.013}));
        return pin;
    };
    auto make_lead = [](const std::string& label, const MAS::Pin& pin, double axial) {
        PinLeadRequest lead;
        lead.label = label;
        lead.pin = pin;
        lead.windowExit = {0.011, axial};
        lead.diameter = 0.0006;
        // Sharp corners: this test is about the run's placement, not its bends (ABT #1296).
        lead.bendRadius = lead.diameter / 2;
        lead.exitX = 0.0;
        return lead;
    };
    SECTION("two strands on one pin too short for both wraps") {
        // 1.5 mm of pin: one wrap needs 0.3 + 2 x 0.6 = 1.5 mm; the second strand's wrap cannot fit.
        auto pin = make_pin("7", 0.0, 0.0015);
        CHECK_THROWS_WITH(OpenMagnetics::Coil::route_leads_to_pins({pin}, {make_lead("strand A", pin, -0.004), make_lead("strand B", pin, -0.005)}, 0.0, 2),
                          Catch::Matchers::ContainsSubstring("Cannot route strand B to pin '7'") &&
                          Catch::Matchers::ContainsSubstring("strand A"));
    }
    SECTION("a pin pitch too narrow to pass between pins") {
        // Pins 1.5 mm apart: a 0.6 mm wire beside pin '2' needs 0.4 + 0.3 + 0.6 = 1.3 mm from both.
        std::vector<MAS::Pin> pins = {make_pin("1", -0.0015, 0.004), make_pin("2", 0.0, 0.004), make_pin("3", 0.0015, 0.004)};
        CHECK_THROWS_WITH(OpenMagnetics::Coil::route_leads_to_pins(pins, {make_lead("the middle lead", pins[1], -0.004)}, 0.0, 2),
                          Catch::Matchers::ContainsSubstring("Cannot route the middle lead to pin '2'") &&
                          Catch::Matchers::ContainsSubstring("pin '1'"));
    }
}

TEST_CASE("With coil_connect_leads_to_pins off (default) autocomplete assigns no pins and keeps the window-exit lead lengths (ABT #1237)",
          "[constructive-model][coil][pins][abt1237][masautocomplete]") {
    SettingsGuard<bool> realWinding(Settings::GetInstance(), &Settings::get_coil_use_real_winding_geometry,
                                    &Settings::set_coil_use_real_winding_geometry, true);
    CHECK_FALSE(Settings::GetInstance().get_coil_connect_leads_to_pins());
    auto build = []() {
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(former_core());
        json coilJson;
        coilJson["bobbin"] = formerBobbinName;
        coilJson["functionalDescription"] = json::parse(R"([
            {"name": "Primary", "numberTurns": 40, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 0.5 - Grade 1"},
            {"name": "Secondary", "numberTurns": 6, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 0.5 - Grade 1"}
        ])");
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
        return magnetic_autocomplete(magnetic);
    };
    auto off = build();
    for (const auto& winding : off.get_coil().get_functional_description()) {
        if (winding.get_connections()) {
            const auto connections = winding.get_connections().value();
            for (const auto& connection : connections) {
                CHECK_FALSE(connection.get_pin_name());
                CHECK_FALSE(connection.get_end());
            }
        }
    }
    auto offCoil = off.get_coil();
    std::vector<ConnectionRoute> offRoutes;
    offCoil.get_connection_reserved_spaces(&offRoutes);
    size_t terminals = 0;
    for (const auto& route : offRoutes) {
        CHECK(route.pinName.empty());
        CHECK(route.pinWaypoints.empty());
        if (route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT) {
            ++terminals;
            // Pre-WP3: the terminal route ends at the window border, its length is its own legs.
            double legs = 0;
            for (size_t k = 0; k + 1 < route.waypoints.size(); ++k) {
                legs += std::hypot(route.waypoints[k + 1][0] - route.waypoints[k][0], route.waypoints[k + 1][1] - route.waypoints[k][1]);
            }
            CHECK_THAT(route.routedLength, Catch::Matchers::WithinAbs(legs, 2e-9));
        }
    }
    CHECK(terminals == 4);
    const auto offLengths = WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(offCoil);

    SettingsGuard<bool> connectLeads(Settings::GetInstance(), &Settings::get_coil_connect_leads_to_pins,
                                     &Settings::set_coil_connect_leads_to_pins, true);
    auto on = build();
    auto onCoil = on.get_coil();
    REQUIRE(onCoil.has_complete_pin_connections());
    const auto onLengths = WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(onCoil);
    // The same winding, longer leads: exactly the runs to the pins.
    REQUIRE(offLengths.size() == onLengths.size());
    for (size_t winding = 0; winding < onLengths.size(); ++winding) {
        CHECK(onLengths[winding][0] > offLengths[winding][0]);
    }
    auto report = check_pin_runs_share_no_copper(onCoil);
    CHECK(report.runs == 4);
}


TEST_CASE("terminal_exit_slots: bundles on consecutive lanes, anchored clear of the plane's routes (ABT #1237)",
          "[constructive-model][coil][pins][abt1237]") {
    // The boost_inductor_complete geometry in miniature (mm -> m): two parallels, entrance rows at the
    // bottom, exit rows at the top, one Z dragback per parallel stepping out at the top of the layer.
    const double od = 0.943e-3;
    auto route = [](const std::string& winding, int64_t parallel, ConnectionKind kind, int side,
                    std::vector<std::vector<double>> waypoints, const std::string& turn) {
        ConnectionRoute r;
        r.winding = winding;
        r.parallel = parallel;
        r.kind = kind;
        r.side = side;
        r.waypoints = waypoints;
        (kind == ConnectionKind::TERMINAL_ENTRANCE ? r.toTurn : r.fromTurn) = turn;
        return r;
    };
    auto bundle = [&](const std::string& winding, int side, bool withDragback) {
        std::vector<ConnectionRoute> routes;
        std::vector<double> attach;
        if (withDragback) {
            routes.push_back(route(winding, 0, ConnectionKind::Z_DRAGBACK, side, {{8.5145e-3, 5.3794e-3}, {9.4575e-3, 5.3794e-3}, {9.4575e-3, -5.3782e-3}}, ""));
            attach.push_back(0);
        }
        // Parallel 0's entrance has no stub, parallel 1's climbs one from its turn: entrance order by span ascending.
        routes.push_back(route(winding, 1, ConnectionKind::TERMINAL_ENTRANCE, side, {{10.795e-3, -6.3231e-3}, {7.5715e-3, -6.3231e-3}, {7.5715e-3, -5.3793e-3}}, "t"));
        attach.push_back(-5.3793e-3);
        routes.push_back(route(winding, 0, ConnectionKind::TERMINAL_ENTRANCE, side, {{10.795e-3, -6.3231e-3}, {7.5715e-3, -6.3231e-3}}, "t"));
        attach.push_back(-6.3231e-3);
        routes.push_back(route(winding, 0, ConnectionKind::TERMINAL_EXIT, side, {{9.4575e-3, 5.3747e-3}, {10.795e-3, 5.3747e-3}}, "t"));
        attach.push_back(5.3747e-3);
        routes.push_back(route(winding, 1, ConnectionKind::TERMINAL_EXIT, side, {{9.4575e-3, 6.3216e-3}, {10.795e-3, 6.3216e-3}}, "t"));
        attach.push_back(6.3216e-3);
        return std::make_pair(routes, attach);
    };
    auto slots_of = [&](const std::vector<ConnectionRoute>& routes, const std::vector<double>& attach) {
        return OpenMagnetics::Coil::terminal_exit_slots(routes, std::vector<double>(routes.size(), od), attach);
    };

    SECTION("entrances at the plane side by side, exits one lane out of their dragback") {
        auto [routes, attach] = bundle("Primary", 0, true);
        auto slots = slots_of(routes, attach);
        CHECK_FALSE(slots[0].has_value());                         // the dragback has no exit slot
        CHECK_THAT(slots[2].value(), Catch::Matchers::WithinAbs(0.0, 1e-15));      // entrance p0
        CHECK_THAT(slots[1].value(), Catch::Matchers::WithinAbs(od, 1e-15));       // entrance p1
        CHECK_THAT(slots[3].value(), Catch::Matchers::WithinAbs(od, 1e-15));       // exit p0
        CHECK_THAT(slots[4].value(), Catch::Matchers::WithinAbs(2 * od, 1e-15));   // exit p1
    }
    SECTION("without the dragback the exits take the plane") {
        auto [routes, attach] = bundle("Primary", 0, false);
        auto slots = slots_of(routes, attach);
        CHECK_THAT(slots[1].value(), Catch::Matchers::WithinAbs(0.0, 1e-15));   // entrance p0
        CHECK_THAT(slots[0].value(), Catch::Matchers::WithinAbs(od, 1e-15));    // entrance p1
        CHECK_THAT(slots[2].value(), Catch::Matchers::WithinAbs(0.0, 1e-15));   // exit p0
        CHECK_THAT(slots[3].value(), Catch::Matchers::WithinAbs(od, 1e-15));    // exit p1
    }
    SECTION("exits order by span descending: a climbing exit takes the lower lane") {
        auto [routes, attach] = bundle("Primary", 0, false);
        routes[3].waypoints = {{9.4575e-3, 7.0e-3}, {9.4575e-3, 8.0e-3}, {10.795e-3, 8.0e-3}};   // parallel 1 exit now climbs
        attach[3] = 7.0e-3;
        auto slots = slots_of(routes, attach);
        CHECK_THAT(slots[3].value(), Catch::Matchers::WithinAbs(0.0, 1e-15));   // exit p1 first
        CHECK_THAT(slots[2].value(), Catch::Matchers::WithinAbs(od, 1e-15));    // exit p0
    }
    SECTION("each isolation side has its own lanes; another winding's link blocks, its own does not") {
        auto [primary, primaryAttach] = bundle("Primary", 0, false);
        auto [secondary, secondaryAttach] = bundle("Secondary", 1, false);
        std::vector<ConnectionRoute> routes = primary;
        std::vector<double> attach = primaryAttach;
        routes.insert(routes.end(), secondary.begin(), secondary.end());
        attach.insert(attach.end(), secondaryAttach.begin(), secondaryAttach.end());
        // A U turnaround of the Primary through the Secondary's exit rows.
        routes.push_back(route("Primary", 0, ConnectionKind::U_ADJACENT, 1, {{9.0e-3, 5.8e-3}, {11.0e-3, 5.8e-3}}, ""));
        attach.push_back(0);
        // ... and the same link on the Primary's own side, through its own exit rows.
        routes.push_back(route("Primary", 0, ConnectionKind::U_ADJACENT, 0, {{9.0e-3, 5.8e-3}, {11.0e-3, 5.8e-3}}, ""));
        attach.push_back(0);
        auto slots = slots_of(routes, attach);
        CHECK_THAT(slots[2].value(), Catch::Matchers::WithinAbs(0.0, 1e-15));        // Primary exit p0 at the plane: own link
        CHECK_THAT(slots[3].value(), Catch::Matchers::WithinAbs(od, 1e-15));
        CHECK_THAT(slots[4 + 1].value(), Catch::Matchers::WithinAbs(0.0, 1e-15));    // Secondary entrance p0: its own lane 0
        CHECK_THAT(slots[4 + 2].value(), Catch::Matchers::WithinAbs(od, 1e-15));     // Secondary exit p0: pushed off the plane
        CHECK_THAT(slots[4 + 3].value(), Catch::Matchers::WithinAbs(2 * od, 1e-15));
    }
    SECTION("inputs that do not match the routes throw") {
        auto [routes, attach] = bundle("Primary", 0, true);
        CHECK_THROWS_WITH(OpenMagnetics::Coil::terminal_exit_slots(routes, {od}, attach), Catch::Matchers::ContainsSubstring("one diameter"));
        attach[1] = std::numeric_limits<double>::quiet_NaN();
        CHECK_THROWS_WITH(slots_of(routes, attach), Catch::Matchers::ContainsSubstring("no attach turn coordinate"));
    }
}

TEST_CASE("route_leads_to_pins starts each run at its exit slot and keeps every run off the other leads' in-window runs (ABT #1237)",
          "[constructive-model][coil][pins][abt1237]") {
    auto make_pin = [](const std::string& name, double x) {
        MAS::Pin pin;
        pin.set_name(name);
        pin.set_shape(PinShape::ROUND);
        pin.set_type(PinDescriptionType::THT);
        pin.set_dimensions({0.0008, 0.0008, 0.004});
        pin.set_coordinates(std::vector<double>({x, -0.014, -0.013}));
        return pin;
    };
    const std::vector<MAS::Pin> pins = {make_pin("A", -0.004), make_pin("B", -0.008)};
    // Lead A leaves at x = 0 lifted 2 mm by the ride-over: its in-window run covers z -11 .. -13 mm.
    // Lead B leaves the same row one OD along, unlifted, at z = -11 mm, and heads -x past x = 0.
    PinLeadRequest a;
    a.label = "lead A";
    a.pin = pins[0];
    a.windowExit = {0.011, -0.004};
    a.diameter = 0.0006;
    // Sharp corners: this test is about the run's placement, not its bends (ABT #1296).
    a.bendRadius = a.diameter / 2;
    a.lift = 0.002;
    a.exitX = 0.0;
    PinLeadRequest b = a;
    b.label = "lead B";
    b.pin = pins[1];
    b.lift = 0.0;
    b.exitX = 0.0006;
    const auto routes = OpenMagnetics::Coil::route_leads_to_pins(pins, {a, b}, 0.0, 2);
    REQUIRE(routes.size() == 2);
    CHECK_THAT(routes[0].waypoints.front()[0], Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(routes[0].waypoints.front()[2], Catch::Matchers::WithinAbs(-0.013, 1e-15));
    CHECK_THAT(routes[1].waypoints.front()[0], Catch::Matchers::WithinAbs(0.0006, 1e-15));
    CHECK_THAT(routes[1].waypoints.front()[2], Catch::Matchers::WithinAbs(-0.011, 1e-15));
    const std::vector<Point3> inWindowA = {{0.0, -0.004, -0.011}, {0.0, -0.004, -0.013}};
    const std::vector<Point3> inWindowB = {{0.0006, -0.004, -0.011}, {0.0006, -0.004, -0.011}};
    CHECK(sampled_distance(routes[1].waypoints, inWindowA) >= 0.0006 - 1e-9);
    CHECK(sampled_distance(routes[0].waypoints, inWindowB) >= 0.0006 - 1e-9);
    CHECK(std::min(sampled_distance(routes[0].waypoints, routes[1].waypoints), sampled_distance(routes[1].waypoints, routes[0].waypoints)) >= 0.0006 - 1e-9);

    SECTION("a lead without an exit slot, or two slots sharing copper, throw") {
        PinLeadRequest noSlot = a;
        noSlot.exitX.reset();
        CHECK_THROWS_WITH(OpenMagnetics::Coil::route_leads_to_pins(pins, {noSlot}, 0.0, 2), Catch::Matchers::ContainsSubstring("has no exit slot x"));
        PinLeadRequest same = b;
        same.exitX = 0.0003;
        CHECK_THROWS_WITH(OpenMagnetics::Coil::route_leads_to_pins(pins, {a, same}, 0.0, 2),
                          Catch::Matchers::ContainsSubstring("The exit slots of lead A and lead B share copper"));
    }
}

// ABT #1172: a pin run's corners are planned for the bend a CONSUMER will draw. MKF plans straight
// legs meeting at right angles around an obstacle edge (the pin rail, here), but MVB++ sweeps each
// corner with a centreline radius R > r, and a rounded corner cuts the INSIDE of the bend — which
// is the edge. With the legs one wire radius off both faces the copper ends up inside the rail:
// measured on the boost PQ 26/25 as 0.001526 mm^3 of overlap. The legs now stand off by
// d >= R - (R - r) sin(theta/2) instead.
TEST_CASE("route_leads_to_pins clears the rail edge for the bend a consumer will draw (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    auto make_pin = [](const std::string& name, double x) {
        MAS::Pin pin;
        pin.set_name(name);
        pin.set_shape(PinShape::ROUND);
        pin.set_type(PinDescriptionType::THT);
        pin.set_dimensions({0.0008, 0.0008, 0.004});
        pin.set_coordinates(std::vector<double>({x, -0.014, -0.013}));
        return pin;
    };
    const std::vector<MAS::Pin> pins = {make_pin("A", -0.004)};
    PinLeadRequest lead;
    lead.label = "lead A";
    lead.pin = pins[0];
    lead.windowExit = {0.011, -0.004};
    lead.diameter = 0.000943;       // the boost PQ 26/25 lead: r = 0.4715 mm
    lead.bendRadius = lead.diameter / 2;
    lead.lift = 0.0;
    lead.exitX = 0.0;
    const double coatedRadius = lead.diameter / 2;

    // A rail the run has to pass under, in the plane of the drop.
    OpenMagnetics::Bobbin::PinRailBlock rail;
    rail.name = "rail 0 block 0";
    rail.centre = {-0.004, -0.0145, -0.0115};
    rail.halfExtents = {0.006, 0.0005, 0.0015};

    SECTION("geometry: the required stand-off is exactly R - (R - r) sin(theta/2)") {
        // Checked against the closed form, independently of the routing code.
        for (double factor : {1.0, 1.05, 1.5}) {
            settings.reset();
            if (factor > 1.0) {
                settings.set_coil_lead_bend_radius_factor(factor);
            }
            const double bendRadius = OpenMagnetics::Settings::resolve_lead_bend_radius(coatedRadius);
            const double clearance = OpenMagnetics::Settings::lead_leg_clearance(coatedRadius, bendRadius, std::numbers::pi / 2);
            CHECK_THAT(bendRadius, Catch::Matchers::WithinRel(factor * coatedRadius, 1e-12));
            CHECK_THAT(clearance,
                       Catch::Matchers::WithinRel(bendRadius - (bendRadius - coatedRadius) / std::sqrt(2.0), 1e-12));
            // The arc of a bend of THIS radius, with legs at THIS stand-off, clears the edge by at
            // least the wire's own radius — which is the whole point.
            const double closestApproach = bendRadius - (bendRadius - clearance) * std::sqrt(2.0);
            CHECK(closestApproach >= coatedRadius - 1e-12);
        }
        settings.reset();
    }

    SECTION("unset: sharp corners, the historical geometry, and the route says so") {
        settings.reset();
        const auto routes = OpenMagnetics::Coil::route_leads_to_pins(pins, {lead}, 0.0, 2, {rail});
        REQUIRE(routes.size() == 1);
        REQUIRE(routes[0].waypoints.size() >= 2);
        // No bend declared -> planned for sharp corners -> the recorded radius is the wire's own,
        // so a consumer that rounds its corners can see MKF did not plan for one.
        CHECK_THAT(routes[0].plannedBendRadius, Catch::Matchers::WithinRel(coatedRadius, 1e-12));
    }

    SECTION("with a bend declared the run moves OUT, and records the radius it was planned for") {
        settings.reset();
        const auto sharpRoutes = OpenMagnetics::Coil::route_leads_to_pins(pins, {lead}, 0.0, 2, {rail});
        REQUIRE(sharpRoutes.size() == 1);
        const double railUnderside = rail.centre[1] - rail.halfExtents[1];

        settings.set_coil_lead_bend_radius_factor(1.05);      // MVB++'s kRoundCornerBendFactor
        PinLeadRequest roundedLead = lead;                    // the radius travels ON the request
        roundedLead.bendRadius = OpenMagnetics::Settings::resolve_lead_bend_radius(coatedRadius);
        const auto roundedRoutes = OpenMagnetics::Coil::route_leads_to_pins(pins, {roundedLead}, 0.0, 2, {rail});
        REQUIRE(roundedRoutes.size() == 1);
        CHECK_THAT(roundedRoutes[0].plannedBendRadius,
                   Catch::Matchers::WithinRel(1.05 * coatedRadius, 1e-12));

        // The planned run MOVES when a bend is declared, and moves by exactly the extra
        // stand-off the arc needs: d(R) - d(r) = (R - r)(1 - 1/sqrt(2)), which is 6.9 um for this
        // lead at R = 1.05 r. Asserted on the waypoints themselves rather than on a distance to a
        // synthetic rail, so the check does not depend on where this test happens to put the rail.
        const double expectedShift = (1.05 * coatedRadius - coatedRadius) * (1 - 1 / std::sqrt(2.0));
        REQUIRE(sharpRoutes[0].waypoints.size() == roundedRoutes[0].waypoints.size());
        double largestShift = 0;
        for (size_t index = 0; index < sharpRoutes[0].waypoints.size(); ++index) {
            for (size_t axis = 0; axis < 3; ++axis) {
                largestShift = std::max(largestShift, std::fabs(roundedRoutes[0].waypoints[index][axis] -
                                                                sharpRoutes[0].waypoints[index][axis]));
            }
        }
        UNSCOPED_INFO("largest waypoint shift " << largestShift * 1e6 << " um, expected "
                      << expectedShift * 1e6 << " um");
        // Reverting the clearance makes this zero: the two runs would be identical.
        CHECK(largestShift > 0);
        CHECK_THAT(largestShift, Catch::Matchers::WithinRel(expectedShift, 1e-9));
        settings.reset();
    }

    SECTION("a bend tighter than the wire is refused, not silently accepted") {
        settings.reset();
        CHECK_THROWS_WITH(settings.set_coil_lead_bend_radius_factor(0.9),
                          Catch::Matchers::ContainsSubstring("cannot be below 1"));
        settings.reset();
    }

    SECTION("a lead planned for a bend tighter than it sweeps is refused (ABT #1296)") {
        settings.reset();
        PinLeadRequest tooTight = lead;
        tooTight.bendRadius = 0.9 * coatedRadius;
        CHECK_THROWS_WITH(OpenMagnetics::Coil::route_leads_to_pins(pins, {tooTight}, 0.0, 2, {rail}),
                          Catch::Matchers::ContainsSubstring("below the radius it sweeps"));
        PinLeadRequest unset = lead;
        unset.bendRadius = 0;                                 // never set: refused, not planned sharp
        CHECK_THROWS_WITH(OpenMagnetics::Coil::route_leads_to_pins(pins, {unset}, 0.0, 2, {rail}),
                          Catch::Matchers::ContainsSubstring("below the radius it sweeps"));
    }
}

TEST_CASE("The bend radius a pin run was planned for reaches the struct a consumer can read (ABT #1172)",
          "[constructive-model][coil][pins][abt1172]") {
    // MVB++ consumes ConnectionRoute from get_connection_layout().routes; PinLeadRoute is
    // route_leads_to_pins' return value and never reaches the enriched magnetic. Settings do not
    // substitute for the field: a consumer that restores its settings after enrichment (or never
    // set them -- the WASM/web path, an older saved design) reads the AMBIENT value at draw time,
    // not what this design was planned with, and would round a corner into a rail whose legs were
    // never offset for it.
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);

    auto inputs = reinforced_offline_inputs();
    auto coil = make_coil({{"Primary", 40, 1, "primary", "Round 0.2 - Grade 1"},
                           {"Secondary", 6, 1, "secondary", "Round 0.2 - Grade 1"}},
                          inputs);
    auto core = former_core();
    REQUIRE_FALSE(coil.assign_pins(coil.resolve_bobbin(), core).skipped);

    const double factor = 1.05;                       // MVB++'s kRoundCornerBendFactor
    settings.set_coil_lead_bend_radius_factor(factor);
    const auto layout = coil.get_connection_layout();
    settings.reset();                                 // as a consumer's SettingsGuard would

    size_t checked = 0;
    for (const auto& route : layout.routes) {
        if (route.pinWaypoints.empty()) {
            continue;                                 // not a run to a pin
        }
        INFO("route " << route.winding << " parallel " << route.parallel << " pin " << route.pinName);
        CHECK(route.plannedBendRadius > 0);
        // The recorded radius is the one the legs were offset for, NOT the wire's own: it must
        // survive the settings being restored, which is the whole reason it lives on the route.
        CHECK(route.plannedBendRadius > 0.5 * factor * 0.0002);
        ++checked;
    }
    REQUIRE(checked > 0);                             // a vacuous pass here would look identical
}

TEST_CASE("A lead is planned for the largest of its wire's, its sleeve's and the drawer's bend radius (ABT #1296)",
          "[constructive-model][coil][pins][abt1296]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto round = find_wire_by_name("Round 0.2 - Grade 1");
    const double coatedRadius = resolve_dimensional_values(round.get_outer_diameter().value()) / 2;
    const double flexibility = WireBend::get_minimum_bend_radius(round, BendCriterion::FLEXIBILITY, BendAxis::ROUND);
    REQUIRE(flexibility > coatedRadius);   // otherwise the next check could not tell the two apart

    SECTION("not a real winding: the pre-#1296 radius, whatever the wire (Alf, 2026-09-21)") {
        CHECK(OpenMagnetics::Coil::lead_bend_radius(round, std::nullopt, coatedRadius, false) == coatedRadius);
        settings.set_coil_lead_bend_radius_factor(1.05);
        CHECK_THAT(OpenMagnetics::Coil::lead_bend_radius(round, std::nullopt, coatedRadius, false),
                   Catch::Matchers::WithinRel(1.05 * coatedRadius, 1e-12));
        settings.reset();
    }

    SECTION("nothing declared: the wire's own IEC 60317-0-1 minimum, not a sharp corner") {
        CHECK_THAT(OpenMagnetics::Coil::lead_bend_radius(round, std::nullopt, coatedRadius, true),
                   Catch::Matchers::WithinRel(flexibility, 1e-12));
    }

    SECTION("a drawer needing more than the wire gets what it declared") {
        const double factor = 1.5 * flexibility / coatedRadius;
        settings.set_coil_lead_bend_radius_factor(factor);
        CHECK_THAT(OpenMagnetics::Coil::lead_bend_radius(round, std::nullopt, coatedRadius, true),
                   Catch::Matchers::WithinRel(factor * coatedRadius, 1e-12));
        settings.reset();
    }

    SECTION("wires no standard covers get buildability alone") {
        OpenMagnetics::Wire litz;
        litz.set_type(WireType::LITZ);
        CHECK(OpenMagnetics::Coil::lead_bend_radius(litz, std::nullopt, 0.0005, true) == 0.0005);
        // IEC 60317-0-1 clause 8.2: no winding test above a 1,600 mm conductor.
        auto thick = find_wire_by_name("Round 2.00 - Grade 1");
        const double thickRadius = resolve_dimensional_values(thick.get_outer_diameter().value()) / 2;
        CHECK(OpenMagnetics::Coil::lead_bend_radius(thick, std::nullopt, thickRadius, true) == thickRadius);
    }

    SECTION("a sleeve whose material rates its bend takes the smallest rated size that holds it") {
        auto material = find_insulation_material_by_name("PTFE extruded tubing");
        auto point = [](double innerDiameter, double wall, double value) {
            MAS::MinimumBendRadiusElement element;
            element.set_inner_diameter(innerDiameter);
            element.set_wall_thickness(wall);
            element.set_value(value);
            return element;
        };
        // Made-up ratings, for the selection rule only. The 1.0 mm size is too narrow for the
        // sleeve, the 1.2 mm / 0.3 mm size is too thin-walled for it, and 1.2 mm / 0.4 mm is the one.
        material.set_minimum_bend_radius(std::vector<MAS::MinimumBendRadiusElement>{
            point(0.0010, 0.0004, 0.008), point(0.0012, 0.0003, 0.009),
            point(0.0012, 0.0004, 0.012), point(0.0020, 0.0004, 0.020)});
        ConnectionSleeve sleeve;
        sleeve.set_material(InsulationMaterialDataOrNameUnion(static_cast<MAS::InsulationMaterial>(material)));
        sleeve.set_inner_diameter(0.0011);
        sleeve.set_wall_thickness(0.0004);
        const double sleeveRadius = 0.0011 / 2 + 0.0004;
        CHECK_THAT(OpenMagnetics::Coil::lead_bend_radius(round, sleeve, sleeveRadius, true),
                   Catch::Matchers::WithinRel(0.012, 1e-12));

        sleeve.set_inner_diameter(0.0025);   // wider than any rated size: refused, never extrapolated
        CHECK_THROWS_WITH(OpenMagnetics::Coil::lead_bend_radius(round, sleeve, 0.0025 / 2 + 0.0004, true),
                          Catch::Matchers::ContainsSubstring("for no tubing"));

        // No rating on the material: the sleeve adds nothing, and the wire rules.
        auto unrated = find_insulation_material_by_name("PTFE extruded tubing");
        REQUIRE_FALSE(unrated.get_minimum_bend_radius());
        sleeve.set_material(InsulationMaterialDataOrNameUnion(static_cast<MAS::InsulationMaterial>(unrated)));
        sleeve.set_inner_diameter(0.0011);
        CHECK_THAT(OpenMagnetics::Coil::lead_bend_radius(round, sleeve, coatedRadius, true),
                   Catch::Matchers::WithinRel(flexibility, 1e-12));
    }
}

TEST_CASE("Every terminal route publishes the bend it was planned for, pin or not (ABT #1296)",
          "[constructive-model][coil][pins][abt1296]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto round = find_wire_by_name("Round 0.2 - Grade 1");
    const double flexibility = WireBend::get_minimum_bend_radius(round, BendCriterion::FLEXIBILITY, BendAxis::ROUND);
    const std::vector<Winding_> windings = {{"Primary", 40, 1, "primary", "Round 0.2 - Grade 1"},
                                            {"Secondary", 6, 1, "secondary", "Round 0.2 - Grade 1"}};

    // Reinforced offline: every lead crosses a margin in a sleeve. No insulation requirements:
    // bare leads, where only the wire's own minimum can move the radius off the sharp corner.
    const bool sleevedCase = GENERATE(true, false);
    INFO((sleevedCase ? "reinforced offline, sleeved leads" : "no insulation inputs, bare leads"));
    auto coil = sleevedCase ? make_coil(windings, reinforced_offline_inputs()) : make_coil(windings);
    auto core = former_core();
    REQUIRE_FALSE(coil.assign_pins(coil.resolve_bobbin(), core).skipped);
    const auto layout = coil.get_connection_layout();

    size_t terminals = 0;
    size_t pinned = 0;
    size_t sleeved = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        INFO("route " << route.winding << " parallel " << route.parallel << " pin '" << route.pinName << "'");
        // Nothing declared by a drawer and no sleeve rating: the wire's own minimum, unless the
        // corner sweeps a sleeve wider than that -- nothing bends tighter than its own radius.
        // A sharp-corner radius on a bare lead is the pre-#1296 behaviour.
        const double swept = route.sleeveOuterDiameter ? route.sleeveOuterDiameter.value() / 2 : 0.0;
        CHECK_THAT(route.plannedBendRadius, Catch::Matchers::WithinRel(std::max(flexibility, swept), 1e-9));
        ++terminals;
        pinned += route.pinWaypoints.empty() ? 0 : 1;
        sleeved += route.sleeveOuterDiameter ? 1 : 0;
    }
    REQUIRE(terminals == 4);
    REQUIRE(pinned > 0);
    // Both branches are really exercised, or the check above could pass on the wrong one.
    CHECK(sleeved == (sleevedCase ? terminals : 0));
    settings.reset();
}

TEST_CASE("An ideal winding's terminal routes keep the corners they always had (ABT #1296)",
          "[constructive-model][coil][pins][abt1296]") {
    // The wire and sleeve bend physics apply ONLY under real winding (Alf, 2026-09-21). An ideal
    // wind with nothing declared still plans sharp corners: the route's radius is the lead's own.
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    REQUIRE_FALSE(settings.get_coil_use_real_winding_geometry());
    auto coil = make_coil({{"Primary", 40, 1, "primary", "Round 0.2 - Grade 1"},
                           {"Secondary", 6, 1, "secondary", "Round 0.2 - Grade 1"}});
    REQUIRE_FALSE(coil.is_real_winding_blocking_applied());
    const auto layout = coil.get_connection_layout();
    auto round = find_wire_by_name("Round 0.2 - Grade 1");
    const double sharp = 0.5 * std::max(round.get_maximum_outer_width(), round.get_maximum_outer_height());
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
        INFO("route " << route.winding << " parallel " << route.parallel);
        CHECK_THAT(route.plannedBendRadius, Catch::Matchers::WithinRel(sharp, 1e-12));
        ++terminals;
    }
    REQUIRE(terminals == 4);
}

TEST_CASE("A stub too short for its two bends is published as a ramp, only under real winding (ABT #1336)",
          "[constructive-model][coil][pins][abt1336]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const std::vector<Winding_> windings = {{"Primary", 40, 2, "primary", "Round 0.2 - Grade 1"},
                                            {"Secondary", 6, 1, "secondary", "Round 0.2 - Grade 1"}};
    // The stub is the axial leg at the TURN end of a terminal route (electrical order: an entrance
    // ends on its turn, an exit starts on it).
    auto stubHeightOf = [](const ConnectionRoute& route) -> std::optional<double> {
        if (route.waypoints.size() < 3) return std::nullopt;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        const auto& turn = entrance ? route.waypoints.back() : route.waypoints.front();
        const auto& next = entrance ? route.waypoints[route.waypoints.size() - 2] : route.waypoints[1];
        const double radial = std::abs(next[0] - turn[0]), axial = std::abs(next[1] - turn[1]);
        if (radial > 1e-12 && axial > 1e-12) return std::nullopt;   // not an axis-aligned stub
        return std::max(radial, axial);
    };

    SECTION("real winding: every short stub is a ramp of length sqrt(h (4R - h)), every long one is not") {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        auto coil = make_coil(windings);
        REQUIRE(coil.is_real_winding_blocking_applied());
        const auto layout = coil.get_connection_layout();
        size_t ramps = 0, straight = 0;
        for (const auto& route : layout.routes) {
            if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
            INFO("route " << route.winding << " parallel " << route.parallel
                          << (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? " entrance" : " exit"));
            const auto h = stubHeightOf(route);
            const double R = route.plannedBendRadius;
            if (route.rampLength) {
                REQUIRE(h);
                CHECK(*h < 2 * R);
                CHECK_THAT(*route.rampLength, Catch::Matchers::WithinRel(std::sqrt(*h * (4 * R - *h)), 1e-12));
                ++ramps;
            }
            else if (h) {
                CHECK(*h >= 2 * R);   // a stub left straight must hold both its bends
                ++straight;
            }
        }
        UNSCOPED_INFO(ramps << " ramps, " << straight << " straight stubs");
        REQUIRE(ramps > 0);   // the second parallel's stub to the shared row is one OD, below 2R
        settings.reset();
    }

    SECTION("an ideal winding publishes no ramp") {
        settings.reset();
        auto coil = make_coil(windings);
        REQUIRE_FALSE(coil.is_real_winding_blocking_applied());
        size_t terminals = 0;
        for (const auto& route : coil.get_connection_layout().routes) {
            if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
            CHECK_FALSE(route.rampLength);
            ++terminals;
        }
        REQUIRE(terminals == 6);
    }
}
