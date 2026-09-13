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
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
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
        // The run ends on the pin base of the assigned pin, and starts at the window exit.
        auto pin = placed_pin(coil, after.pinName);
        const auto& pinEnd = after.kind == ConnectionKind::TERMINAL_ENTRANCE ? after.pinWaypoints.front() : after.pinWaypoints.back();
        CHECK_THAT(pinEnd[0], Catch::Matchers::WithinAbs(pin.base[0], 1e-12));
        CHECK_THAT(pinEnd[1], Catch::Matchers::WithinAbs(pin.base[1], 1e-12));
        CHECK_THAT(pinEnd[2], Catch::Matchers::WithinAbs(pin.base[2], 1e-12));
        const auto& windowEnd = after.kind == ConnectionKind::TERMINAL_ENTRANCE ? after.pinWaypoints.back() : after.pinWaypoints.front();
        const auto& borderEnd = after.kind == ConnectionKind::TERMINAL_ENTRANCE ? after.waypoints.front() : after.waypoints.back();
        // The border point (radial, axial) sits on the front face at (0, axial, -radial).
        CHECK_THAT(windowEnd[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
        CHECK_THAT(windowEnd[1], Catch::Matchers::WithinAbs(borderEnd[1], 1e-12));
        CHECK_THAT(windowEnd[2], Catch::Matchers::WithinAbs(-borderEnd[0], 1e-12));
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

TEST_CASE("route_lead_to_pin walks the flange face to the pin base (ABT #1172)", "[constructive-model][coil][pins][abt1172]") {
    MAS::Pin vertical;
    vertical.set_name("1");
    vertical.set_shape(PinShape::ROUND);
    vertical.set_type(PinDescriptionType::THT);
    vertical.set_dimensions({0.001, 0.001, 0.004});
    vertical.set_coordinates(std::vector<double>({-0.005, -0.012, 0.006}));
    auto route = OpenMagnetics::Coil::route_lead_to_pin(vertical, {0.010, 0.003}, 0.0);
    REQUIRE(route.waypoints.size() == 4);
    // Exit on the front face (0, 0.003, -0.010); down to the flange face y = -0.010 (0.013), across
    // to the row z = 0.006 (0.016), along the row to x = -0.005 (0.005).
    CHECK_THAT(route.length, Catch::Matchers::WithinAbs(0.013 + 0.016 + 0.005, 1e-12));
    CHECK_THAT(route.waypoints.front()[2], Catch::Matchers::WithinAbs(-0.010, 1e-12));
    CHECK_THAT(route.waypoints.back()[0], Catch::Matchers::WithinAbs(-0.005, 1e-12));
    CHECK_THAT(route.waypoints.back()[1], Catch::Matchers::WithinAbs(-0.010, 1e-12));
    CHECK_THAT(route.waypoints.back()[2], Catch::Matchers::WithinAbs(0.006, 1e-12));

    MAS::Pin horizontal = vertical;
    horizontal.set_coordinates(std::vector<double>({0.004, 0.010, -0.020}));
    horizontal.set_rotation(std::vector<double>({90, 0, 0}));
    route = OpenMagnetics::Coil::route_lead_to_pin(horizontal, {0.010, 0.003}, 0.0);
    // On along -Z to the base plane z = -0.018 (0.008), up to the row y = 0.010 (0.007), along the
    // row to x = 0.004 (0.004).
    CHECK_THAT(route.length, Catch::Matchers::WithinAbs(0.008 + 0.007 + 0.004, 1e-12));
}

TEST_CASE("Autocomplete assigns pins when the bobbin has them (ABT #1172)", "[constructive-model][coil][pins][abt1172][masautocomplete]") {
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

