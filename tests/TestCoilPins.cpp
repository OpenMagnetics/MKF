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


// ABT #1336 (Alf, 2026-10-19): a terminal stub's height where its lead REALLY leaves its turn --
// pitch-true, at the lead's exit slot -- not at the turn's 2-D station. On a round column the turn
// passes slot x at the azimuth asin(x / r) off the connection plane and has advanced by its own
// per-revolution advance (to the next turn for an entrance, from the previous for an exit) over
// that angle. nullopt for a route with no axis-aligned axial stub.
struct PitchTrueStub {
    double station = 0;    // |row - station|: the 2-D height
    double pitchTrue = 0;  // |row - where the turn is at the slot|
};
std::optional<PitchTrueStub> pitchTrueStubOf(const OpenMagnetics::Coil& coil, const ConnectionRoute& route) {
    if (route.waypoints.size() < 3) return std::nullopt;
    const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
    const auto& turnEnd = entrance ? route.waypoints.back() : route.waypoints.front();
    const auto& next = entrance ? route.waypoints[route.waypoints.size() - 2] : route.waypoints[1];
    if (std::abs(next[0] - turnEnd[0]) > 1e-12 || !(std::abs(next[1] - turnEnd[1]) > 1e-12)) return std::nullopt;
    REQUIRE(route.exitSlot);
    const auto turns = coil.get_turns_description().value();
    const std::string& turnName = entrance ? route.toTurn : route.fromTurn;
    size_t at = turns.size();
    for (size_t k = 0; k < turns.size(); ++k)
        if (turns[k].get_name() == turnName) at = k;
    REQUIRE(at < turns.size());
    double advance = 0;
    const auto& turn = turns[at];
    if (entrance) {
        for (size_t k = at + 1; k < turns.size(); ++k)
            if (turns[k].get_winding() == turn.get_winding() && turns[k].get_parallel() == turn.get_parallel()) {
                advance = turns[k].get_coordinates()[1] - turn.get_coordinates()[1];
                break;
            }
    }
    else {
        for (size_t k = at; k-- > 0;)
            if (turns[k].get_winding() == turn.get_winding() && turns[k].get_parallel() == turn.get_parallel()) {
                advance = turn.get_coordinates()[1] - turns[k].get_coordinates()[1];
                break;
            }
    }
    const double x = route.exitSlot.value();
    REQUIRE(std::abs(x) <= turnEnd[0]);
    const double turnY = turnEnd[1] - advance * std::asin(x / turnEnd[0]) / (2 * std::numbers::pi);
    return PitchTrueStub{std::abs(next[1] - turnEnd[1]), std::abs(next[1] - turnY)};
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
// ABT #1237: the boost inductor on PQ 26/25 pins (ABT #2268: its pin runs are parked; the assignment stays)
// -----------------------------------------------------------------------------------------

namespace {

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

TEST_CASE("The boost inductor's four ends on PQ 26/25 pins: strands on adjacent inner pins, starts towards -X (ABT #1237)",
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
    // Strands on adjacent inner pins, starts towards -X. ABT #2268: the runs of these leads to their
    // pins are parked; only the assignment written into the connections is checked.
    CHECK(connection_of(coil, "Primary", End::START, 0).get_pin_name().value() == "2");
    CHECK(connection_of(coil, "Primary", End::START, 1).get_pin_name().value() == "3");
    CHECK(connection_of(coil, "Primary", End::FINISH, 0).get_pin_name().value() == "4");
    CHECK(connection_of(coil, "Primary", End::FINISH, 1).get_pin_name().value() == "5");
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

    SettingsGuard<bool> connectLeads(Settings::GetInstance(), &Settings::get_coil_connect_leads_to_pins,
                                     &Settings::set_coil_connect_leads_to_pins, true);
    auto on = build();
    auto onCoil = on.get_coil();
    REQUIRE(onCoil.has_complete_pin_connections());
    // ABT #2268: the runs of the leads to their pins are parked; with the flag on only the assignment
    // reaches the design.
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
        return OpenMagnetics::Coil::terminal_exit_slots(routes, std::vector<double>(routes.size(), od), attach, false,
                                                        std::vector<double>(routes.size(), 0.0), std::vector<std::optional<TerminalTurnPath>>{}, {});
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
        CHECK_THROWS_WITH(OpenMagnetics::Coil::terminal_exit_slots(routes, {od}, attach, false, std::vector<double>(routes.size(), 0.0), {}, {}), Catch::Matchers::ContainsSubstring("one diameter"));
        attach[1] = std::numeric_limits<double>::quiet_NaN();
        CHECK_THROWS_WITH(slots_of(routes, attach), Catch::Matchers::ContainsSubstring("no attach turn coordinate"));
    }
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
        sleeved += route.sleeveOuterDiameter ? 1 : 0;
    }
    REQUIRE(terminals == 4);
    // Both branches are really exercised, or the check above could pass on the wrong one.
    CHECK(sleeved == (sleevedCase ? terminals : 0));
    settings.reset();
}

namespace {
// Spec section 12: the face a pin's lead leaves on (its rail's side of the column) and the half of that
// face its pin's x lies on, in the face's own frame (face 1 is the bobbin frame turned 180 degrees).
int pin_face(const MAS::Pin& pin) {
    REQUIRE(pin.get_coordinates());
    REQUIRE(std::abs(pin.get_coordinates().value()[2]) > 1e-9);
    return pin.get_coordinates().value()[2] < 0 ? 0 : 1;
}
int pin_half(const MAS::Pin& pin) {
    const double x = pin.get_coordinates().value()[0];
    const double faceX = pin_face(pin) == 0 ? x : -x;
    REQUIRE(std::abs(faceX) > 1e-9);
    return faceX < 0 ? -1 : 1;
}
// ABT #2268: the pin a terminal route's end is assigned to, read from the coil's connections (the
// assignment in MAS). The route itself names no pin while the runs to the pins are parked.
std::string assigned_pin_of(const OpenMagnetics::Coil& coil, const ConnectionRoute& route) {
    const End end = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? End::START : End::FINISH;
    return connection_of(coil, route.winding, end).get_pin_name().value();
}
}  // namespace

TEST_CASE("Pins decide the side: each pinned winding leaves on the face toward its own rail, on its pin's half (spec section 12, ABT #1336)",
          "[constructive-model][coil][pins][abt1336][abt1172]") {
    // The ABT #1172 pinned flyback, built exactly as "The bend radius a pin run was planned for reaches
    // the struct a consumer can read" built it (that test was removed with the parked pin runs, ABT #2268):
    // the Primary's pins on rail 0 (z < 0), the Secondary's on rail 1.
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
    settings.reset();
    auto bobbin = coil.resolve_bobbin();

    std::map<std::string, int> faceOfWinding;
    std::set<int> faces;
    size_t pinnedTerminals = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        const auto pinName = assigned_pin_of(coil, route);
        INFO("route " << route.winding << " parallel " << route.parallel << " pin '" << pinName << "'");
        REQUIRE(route.exitSlot.has_value());
        const auto pin = bobbin.get_pin(pinName);
        const int face = pin_face(pin);
        // (a) the lead leaves on the face toward its pin's rail.
        CHECK(route.side == face);
        INFO("exit slot " << route.exitSlot.value() * 1e3 << " mm, half " << route.laneHalf.value_or(0));
        // (b) the lead packs on the half of its face its pin's x lies on.
        REQUIRE(route.laneHalf.has_value());
        CHECK(route.laneHalf.value() == pin_half(pin));
        CHECK(route.exitSlot.value() * route.laneHalf.value() >= -1e-12);
        faceOfWinding[route.winding] = face;
        faces.insert(face);
        ++pinnedTerminals;
    }
    REQUIRE(pinnedTerminals == 4);
    // Both faces are exercised, or the face check could pass on the one-face default.
    REQUIRE(faces.size() == 2);
    // Every route of a pinned winding, its dragbacks too, is on its winding's face.
    for (const auto& route : layout.routes) {
        INFO("route " << route.winding << " parallel " << route.parallel << " kind " << int(route.kind));
        REQUIRE(faceOfWinding.count(route.winding));
        CHECK(route.side == faceOfWinding.at(route.winding));
    }
}

TEST_CASE("Pins decide the side with coil_connect_leads_to_pins off too: the same faces, halves and lanes (spec section 12, ABT #1336)",
          "[constructive-model][coil][pins][abt1336][abt1237][masautocomplete]") {
    // The fixture of "With coil_connect_leads_to_pins off (default) autocomplete assigns no pins ...", exactly.
    SettingsGuard<bool> realWinding(Settings::GetInstance(), &Settings::get_coil_use_real_winding_geometry,
                                    &Settings::set_coil_use_real_winding_geometry, true);
    REQUIRE_FALSE(Settings::GetInstance().get_coil_connect_leads_to_pins());
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
    auto offCoil = build().get_coil();
    REQUIRE_FALSE(offCoil.has_complete_pin_connections());
    std::vector<ConnectionRoute> offRoutes;
    offCoil.get_connection_reserved_spaces(&offRoutes);

    SettingsGuard<bool> connectLeads(Settings::GetInstance(), &Settings::get_coil_connect_leads_to_pins,
                                     &Settings::set_coil_connect_leads_to_pins, true);
    auto onCoil = build().get_coil();
    REQUIRE(onCoil.has_complete_pin_connections());
    std::vector<ConnectionRoute> onRoutes;
    onCoil.get_connection_reserved_spaces(&onRoutes);

    REQUIRE(offRoutes.size() == onRoutes.size());
    std::set<int> faces;
    size_t terminals = 0;
    for (size_t index = 0; index < onRoutes.size(); ++index) {
        const auto& on = onRoutes[index];
        const auto& off = offRoutes[index];
        INFO("route " << on.winding << " parallel " << on.parallel << " kind " << int(on.kind));
        REQUIRE(off.winding == on.winding);
        REQUIRE(off.parallel == on.parallel);
        REQUIRE(off.kind == on.kind);
        CHECK(off.side == on.side);
        CHECK(off.laneHalf == on.laneHalf);
        CHECK(off.laneBundleRank == on.laneBundleRank);
        REQUIRE(off.exitSlot.has_value() == on.exitSlot.has_value());
        if (on.exitSlot) {
            CHECK_THAT(off.exitSlot.value(), Catch::Matchers::WithinAbs(on.exitSlot.value(), 1e-12));
        }
        if (on.kind == ConnectionKind::TERMINAL_ENTRANCE || on.kind == ConnectionKind::TERMINAL_EXIT) {
            ++terminals;
            CHECK(off.pinWaypoints.empty());            // off: the lead ends in free air
            // ABT #2268: on, too: the runs to the pins are parked, so the assigned lead ends in free air.
            CHECK(on.pinWaypoints.empty());
            // The face is the pin's, flag or not.
            CHECK(on.side == pin_face(onCoil.resolve_bobbin().get_pin(assigned_pin_of(onCoil, on))));
            faces.insert(on.side);
        }
    }
    CHECK(terminals == 4);
    REQUIRE(faces.size() == 2);   // the pins put the windings on both faces: the comparison is not the one-face default
}

TEST_CASE("A lead sleeve ends where the layer stack starts, and the terminal space is the bare wire (spec section 15, ABT #1336)",
          "[constructive-model][coil][pins][abt1336][abt1172]") {
    // Alf, 2026-10-26: "The sleeve goes only until the winding layers start; it is designed to maintain
    // insulation around the bobbin and the core, not between the wires, so the terminal space should be
    // reserved for the normal wire, without sleeve." The ABT #1172 pinned flyback, built exactly as "The
    // bend radius a pin run was planned for reaches the struct a consumer can read" (reinforced offline
    // inputs: every lead is sleeved, 0.32 mm bore + 0.41 mm wall = 1.14 mm, 0.22 mm into the winding).
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
    settings.reset();

    const auto wire = find_wire_by_name("Round 0.2 - Grade 1");
    const double bare = std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height());

    // (a) The sleeve's planned span (Alf, 2026-10-26, "where the layer stack starts"): from the window
    // border along the lead to where it first meets a wound layer, conduction or insulation, plus the
    // sleeve's own overlapIntoWinding -- walked here independently of the plan, on the coil's layers.
    std::vector<std::array<double, 4>> stack;
    const auto layers = coil.get_layers_description().value();   // by value (the ABT #650 class)
    for (const auto& layer : layers) {
        const auto c = layer.get_coordinates();
        const auto d = layer.get_dimensions();
        REQUIRE(c.size() >= 2);
        REQUIRE(d.size() >= 2);
        stack.push_back({c[0] - d[0] / 2, c[0] + d[0] / 2, c[1] - d[1] / 2, c[1] + d[1] / 2});
    }
    auto inStack = [&](const std::vector<double>& p) {
        for (const auto& r : stack) {
            if (p[0] > r[0] + 1e-12 && p[0] < r[1] - 1e-12 && p[1] > r[2] + 1e-12 && p[1] < r[3] - 1e-12) {
                return true;
            }
        }
        return false;
    };
    std::optional<double> lowFace, highFace;
    std::optional<double> sleeveDiameter;
    std::map<std::string, size_t> entered;            // winding -> leads that run in the margin band first
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        INFO("route " << route.winding << (entrance ? " entrance" : " exit"));
        REQUIRE(route.sleeveOuterDiameter.has_value());
        REQUIRE(route.marginInnerFace.has_value());
        REQUIRE(route.sleeveWindingEntry.has_value());
        REQUIRE(route.sleeveSpan.has_value());
        sleeveDiameter = route.sleeveOuterDiameter.value();
        const double face = route.marginInnerFace.value();
        (face > 0 ? highFace : lowFace) = face;
        std::optional<double> overlap;
        // By value: get_functional_description() and get_connections() return by value (the ABT #650 class).
        const auto winding = coil.get_functional_description()[coil.get_winding_index_by_name(route.winding)];
        const auto connections = winding.get_connections().value();
        for (const auto& connection : connections) {
            if (connection.get_end() && connection.get_end().value() == (entrance ? End::START : End::FINISH) && connection.get_sleeve()) {
                overlap = connection.get_sleeve()->get_overlap_into_winding();
            }
        }
        REQUIRE(overlap.has_value());
        auto path = route.waypoints;
        if (!entrance) {
            std::reverse(path.begin(), path.end());
        }
        // March along the lead in 0.1 um steps until a point lies inside a layer.
        double entry = 0;
        bool met = false;
        for (size_t k = 0; k + 1 < path.size() && !met; ++k) {
            const double length = std::hypot(path[k + 1][0] - path[k][0], path[k + 1][1] - path[k][1]);
            const size_t steps = size_t(std::ceil(length / 1e-7));
            for (size_t m = 1; m <= steps && !met; ++m) {
                const double t = double(m) / double(steps);
                if (inStack({path[k][0] + t * (path[k + 1][0] - path[k][0]), path[k][1] + t * (path[k + 1][1] - path[k][1])})) {
                    entry += t * length;
                    met = true;
                }
            }
            if (!met) {
                entry += length;
            }
        }
        REQUIRE(met);
        if (entry > 0) {
            ++entered[route.winding];
        }
        INFO("entry " << entry * 1e3 << " mm, overlap " << overlap.value() * 1e3 << " mm, planned entry "
                      << route.sleeveWindingEntry.value() * 1e3 << " mm, span " << route.sleeveSpan.value() * 1e3 << " mm");
        // The march lands within one 0.1 um step past the exact face (it takes the first point strictly inside).
        CHECK(entry - route.sleeveWindingEntry.value() >= -1e-12);
        CHECK(entry - route.sleeveWindingEntry.value() <= 1e-7 + 1e-12);
        CHECK_THAT(route.sleeveSpan.value(), Catch::Matchers::WithinAbs(route.sleeveWindingEntry.value() + overlap.value(), 1e-12));
        ++terminals;
    }
    REQUIRE(terminals == 4);
    REQUIRE(lowFace.has_value());
    REQUIRE(highFace.has_value());
    REQUIRE(sleeveDiameter.has_value());
    REQUIRE(sleeveDiameter.value() > bare);           // or the sizes could not be told apart
    // Both windings' leads run in from the border through empty window before they meet the stack (the
    // Primary along its margin-band row, the Secondary at its own height), for every winding alike.
    CHECK(entered["Primary"] == 2);
    CHECK(entered["Secondary"] == 2);

    // (b) Every terminal reservation is the bare wire across its run: the margin-band edge run lies beside
    // the tape layers, where there is no sleeve.
    size_t inside = 0, inBandOnly = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.isTerminal || !space.layer.empty()) {
            continue;
        }
        const double low = space.coordinates[1] - space.dimensions[1] / 2;
        const double high = space.coordinates[1] + space.dimensions[1] / 2;
        const double across = std::min(space.dimensions[0], space.dimensions[1]);
        INFO("space " << space.winding << " kind " << int(space.kind) << " y [" << low * 1e3 << ", " << high * 1e3 << "] mm, across "
                      << across * 1e3 << " mm");
        CHECK_THAT(across, Catch::Matchers::WithinAbs(bare, 1e-12));
        ++(high > lowFace.value() + 1e-9 && low < highFace.value() - 1e-9 ? inside : inBandOnly);
    }
    CHECK(inside == 4);       // the Primary's two stubs, the Secondary's two radial runs
    CHECK(inBandOnly == 2);   // the Primary's two margin-band edge runs, bare too
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

    SECTION("real winding: every short stub is a cosine-S ramp whose curvature stays within the bend, every long one is not") {
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
            const auto stub = pitchTrueStubOf(coil, route);
            const auto h = stub ? std::optional<double>(stub->pitchTrue) : std::optional<double>{};
            const double R = route.plannedBendRadius;
            if (route.rampLength) {
                REQUIRE(h);
                CHECK(*h < 2 * R);
                // The cosine S at the in-surface radius that leaves room for the column's own
                // curvature: former_core is a PQ, a ROUND column, so rho is the turn's radius.
                const auto& turnEnd = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.waypoints.back() : route.waypoints.front();
                const double rho = turnEnd[0];
                const double Rg = 1.0 / std::sqrt(1.0 / (R * R) - 1.0 / (rho * rho));
                CHECK_THAT(*route.rampLength, Catch::Matchers::WithinRel(std::numbers::pi * std::sqrt(*h * Rg / 2), 1e-12));
                // ... and the S's peak curvature is then exactly 1/Rg, never tighter than the bend
                CHECK_THAT(*h * std::numbers::pi * std::numbers::pi / (2 * *route.rampLength * *route.rampLength),
                           Catch::Matchers::WithinRel(1.0 / Rg, 1e-12));
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

TEST_CASE("An edge-wound rectangular secondary held whole by its section leaves by local radial leads at its own turns' height, with no ramp (ABT #1481, #1487)",
          "[constructive-model][coil][pins][abt1336][abt1481]") {
    // The web's PQ 65/60 with a 10 x 0.8 mm rectangular secondary lying with its width radial, so
    // wound EDGEWISE round the column: its edgewise minimum (about 25 mm) exceeds the turn's own
    // radius. This layout once blanked the web's real-winding views (ABT #1481: the leads' ramps
    // were checked at the edgewise radius). Under ABT #1487 R1 the secondary's section holds each
    // of its parallels whole, so every one of its terminals is LOCAL to that section; with one
    // layer there is nothing of its own section outward of the turn to cross, so each lead leaves
    // radially at its own turn's height: no stub, and so no ramp. The ramp rule of #1481 itself is
    // tested directly, on the strip's own numbers, in the next test case.
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    const json rectangularJson = json::parse(R"({"type": "rectangular", "standard": "IEC 60317", "material": "copper",
        "numberConductors": 1, "conductingWidth": {"nominal": 0.01}, "conductingHeight": {"nominal": 0.0008},
        "outerWidth": {"nominal": 0.010085}, "outerHeight": {"nominal": 0.000885},
        "coating": {"type": "enamelled", "grade": 1, "breakdownVoltage": 750}})");
    OpenMagnetics::Wire rectangular(rectangularJson);
    // As the web lays it out: sections stacked axially, turns concentric inside each.
    auto coil = OpenMagneticsTesting::get_quick_coil({45, 3}, {1, 4}, "PQ 65/60", 1,
                                                     WindingOrientation::CONTIGUOUS, WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED, CoilAlignment::CENTERED,
                                                     {find_wire_by_name("Round 0.5 - Grade 1"), rectangular});
    // Wound by consecutive turns, as the web does, so each parallel's end turn is off its own row.
    const std::string secondaryName = coil.get_functional_description()[1].get_name();
    coil.preload_winding_style_overrides({{secondaryName, WindingStyle::WIND_BY_CONSECUTIVE_TURNS}});
    REQUIRE(coil.wind());
    REQUIRE(coil.is_real_winding_blocking_applied());

    const double flatwise = WireBend::get_flexibility_bend_radius_if_standardised(rectangular, BendAxis::FLATWISE).value();
    const double edgewise = WireBend::get_flexibility_bend_radius_if_standardised(rectangular, BendAxis::EDGEWISE).value();
    REQUIRE(flatwise < edgewise);
    const double coatedRadius = 0.5 * std::max(rectangular.get_maximum_outer_width(), rectangular.get_maximum_outer_height());

    // R1's premise: each secondary parallel's turns all lie in one section.
    const auto turns = coil.get_turns_description().value();
    std::map<int64_t, std::set<std::string>> sectionsByParallel;
    for (const auto& turn : turns) {
        if (turn.get_winding() == secondaryName) {
            sectionsByParallel[turn.get_parallel()].insert(turn.get_section().value());
        }
    }
    REQUIRE(sectionsByParallel.size() == 4);
    for (const auto& [parallel, sections] : sectionsByParallel) {
        INFO("Secondary parallel " << parallel);
        REQUIRE(sections.size() == 1);
    }

    OpenMagnetics::ConnectionLayout layout;
    REQUIRE_NOTHROW(layout = coil.get_connection_layout());
    std::set<std::pair<int64_t, bool>> terminals;   // (parallel, entrance)
    for (const auto& route : layout.routes) {
        if (route.winding != secondaryName) continue;
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        INFO("Secondary parallel " << route.parallel << (entrance ? " entrance" : " exit"));
        terminals.insert({route.parallel, entrance});
        // Local and radial: the turn and the window exit, nothing between -- no stub, so no ramp.
        REQUIRE(route.waypoints.size() == 2);
        CHECK_FALSE(route.rampLength.has_value());
        const auto& turnEnd = entrance ? route.waypoints.back() : route.waypoints.front();
        const auto& windowEnd = entrance ? route.waypoints.front() : route.waypoints.back();
        const std::string turnName = entrance ? route.toTurn : route.fromTurn;
        const auto turn = std::find_if(turns.begin(), turns.end(), [&](const Turn& t) { return t.get_name() == turnName; });
        REQUIRE(turn != turns.end());
        CHECK(turn->get_winding() == secondaryName);
        CHECK(turn->get_parallel() == route.parallel);
        CHECK(turn->get_section().value() == *sectionsByParallel.at(route.parallel).begin());
        // The lead leaves at its own turn's height, straight outward.
        CHECK(turnEnd[0] == turn->get_coordinates()[0]);
        CHECK(turnEnd[1] == turn->get_coordinates()[1]);
        CHECK(windowEnd[1] == turnEnd[1]);
        CHECK(windowEnd[0] > turnEnd[0]);
        // Its corners are planned for its FLATWISE limit and its coated radius only, below the
        // edgewise one (ABT #2012, Alf 2026-10-23: the IEC 60317-0-2 Table 6 edgewise mandrel qualifies
        // the coating, it is not a forming limit; analyze_results reports the lead instead).
        CHECK(route.plannedBendRadius < edgewise);
        CHECK_THAT(route.plannedBendRadius, Catch::Matchers::WithinRel(std::max(coatedRadius, flatwise), 1e-12));
    }
    REQUIRE(terminals.size() == 8);   // an entrance and an exit for each of the 4 parallels
    settings.reset();
}

TEST_CASE("A rectangular lead's ramp is planned at the in-surface axis's minimum, not the tighter of both (ABT #1481)",
          "[constructive-model][coil][pins][abt1481]") {
    // The S of a ramp lies ON the turn's surface and bends the wire about the RADIAL axis; the
    // column's curvature is the turn's own bend, already wound. So a rectangular lead's ramp is
    // checked at Rg = max(isotropic coated radius with the column curvature taken out, the
    // flatwise minimum), never at the edgewise minimum. Fed the 10 x 0.8 mm strip's own numbers.
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    const json rectangularJson = json::parse(R"({"type": "rectangular", "standard": "IEC 60317", "material": "copper",
        "numberConductors": 1, "conductingWidth": {"nominal": 0.01}, "conductingHeight": {"nominal": 0.0008},
        "outerWidth": {"nominal": 0.010085}, "outerHeight": {"nominal": 0.000885},
        "coating": {"type": "enamelled", "grade": 1, "breakdownVoltage": 750}})");
    OpenMagnetics::Wire rectangular(rectangularJson);
    const double flatwise = WireBend::get_flexibility_bend_radius_if_standardised(rectangular, BendAxis::FLATWISE).value();
    const double edgewise = WireBend::get_flexibility_bend_radius_if_standardised(rectangular, BendAxis::EDGEWISE).value();
    REQUIRE(flatwise < edgewise);
    // The coated half-size c, the radius the lead's corner sweeps; nothing declared, no sleeve, so
    // the isotropic part of the lead's limit is c itself.
    const double c = 0.5 * std::max(rectangular.get_maximum_outer_width(), rectangular.get_maximum_outer_height());
    const double isotropicRadius = OpenMagnetics::Coil::lead_isotropic_bend_radius(std::nullopt, c, true);
    REQUIRE(isotropicRadius == c);

    for (const double rho : {0.006, 0.010, 0.0125, 0.020, 0.050}) {
        INFO("turn radius " << rho * 1e3 << " mm");
        const double Rg = OpenMagnetics::Coil::lead_ramp_governing_radius(rectangular, isotropicRadius, 1.0 / rho);
        const double isotropic = 1.0 / std::sqrt(1.0 / (c * c) - 1.0 / (rho * rho));
        CHECK_THAT(Rg, Catch::Matchers::WithinRel(std::max(isotropic, flatwise), 1e-12));
        CHECK(Rg < edgewise);
        for (const double h : {0.0005, 0.002, 0.01}) {
            INFO("stub height " << h * 1e3 << " mm");
            CHECK_THAT(OpenMagnetics::Coil::lead_ramp_length(h, Rg), Catch::Matchers::WithinRel(std::numbers::pi * std::sqrt(h * Rg / 2), 1e-12));
        }
    }
    // A flat face adds no curvature: Rg is the larger of c and the flatwise minimum.
    CHECK_THAT(OpenMagnetics::Coil::lead_ramp_governing_radius(rectangular, isotropicRadius, 0.0),
               Catch::Matchers::WithinRel(std::max(c, flatwise), 1e-12));
    // A column curved tighter than the lead may bend leaves no ramp that can climb.
    CHECK_THROWS(OpenMagnetics::Coil::lead_ramp_governing_radius(rectangular, isotropicRadius, 1.0 / (0.9 * c)));
    settings.reset();
}

TEST_CASE("Pinned sibling leads stand far enough apart for the ramps between them (ABT #1336)",
          "[constructive-model][coil][pins][abt1336]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto coil = make_coil({{"Primary", 40, 2, "primary", "Round 0.2 - Grade 1"},
                           {"Secondary", 6, 1, "secondary", "Round 0.2 - Grade 1"}});
    auto core = former_core();
    REQUIRE_FALSE(coil.assign_pins(coil.resolve_bobbin(), core).skipped);
    const auto layout = coil.get_connection_layout();

    // The window exit of each pinned terminal, at its exit slot in the bobbin frame (face 1 is the bobbin
    // frame turned 180 degrees about the column axis). ABT #2268: read from the route's exit slot, where
    // the parked pin run used to start (its first point's x was exactly this, spec section 12).
    struct Lead { double x; double span; };
    auto wire = find_wire_by_name("Round 0.2 - Grade 1");
    const double coatedDiameter = std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height());
    std::map<std::pair<std::string, bool>, std::vector<Lead>> groups;   // (winding, entrance) -> by parallel
    size_t ramped = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        REQUIRE(route.exitSlot.has_value());
        REQUIRE((route.side == 0 || route.side == 1));
        const double exitX = route.side == 0 ? route.exitSlot.value() : -route.exitSlot.value();
        const double span = route.rampLength ? *route.rampLength + route.plannedBendRadius : 0.0;
        ramped += route.rampLength ? 1 : 0;
        auto& group = groups[{route.winding, entrance}];
        if (group.size() <= size_t(route.parallel)) group.resize(route.parallel + 1);
        group[route.parallel] = {exitX, span};
    }
    REQUIRE(ramped > 0);   // otherwise this checks the plain lane grid only
    size_t pairs = 0;
    for (const auto& [key, leads] : groups) {
        for (size_t k = 1; k < leads.size(); ++k) {
            INFO(key.first << (key.second ? " entrance" : " exit") << " parallels " << k - 1 << " and " << k);
            // The ramp between two siblings (lane_footprint): an entrance's ramp runs from its lane
            // towards lower x, an exit's towards higher x, so the one between them is the higher
            // sibling's for an entrance and the lower sibling's for an exit. ABT #1705 (re-derived
            // with Alf's approval of 2026-10-24): entrances now pack outwards on x <= 0, so the
            // pair is ordered by x, not by parallel index.
            const auto& lower = leads[k].x < leads[k - 1].x ? leads[k] : leads[k - 1];
            const auto& higher = leads[k].x < leads[k - 1].x ? leads[k - 1] : leads[k];
            const double between = key.second ? higher.span : lower.span;
            const double gap = std::abs(leads[k].x - leads[k - 1].x);
            CHECK(gap > between);   // the ramp fits between them ...
            CHECK(gap - between >= coatedDiameter - 1e-12);   // ... with a coated wire's centreline pitch to spare
            ++pairs;
        }
    }
    REQUIRE(pairs > 0);
    settings.reset();
}

// ABT #1336 (Alf, 2026-10-19: "MKF plans pitch-true"). 14_dab's Primary parallel 1 EXIT leaves a final
// landing that descends a whole window per revolution: at its 2-D station the stub to the edge row is
// 0.934 mm, below the 2R its two bends need, so the station rule planned a ramp -- but at its exit slot
// (8.66 mm) the turn is 3.77 mm from the row, which holds both bends, and the drawn lead attaches there.
// Every stub is planned from where its own turn is at its own slot: short ones ramp at the cosine S of
// that height, long ones stay straight.
TEST_CASE("A terminal stub ramps or not by its height at its exit slot, not at its station (ABT #1336)",
          "[constructive-model][coil][pins][abt1336]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    clear_databases();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);   // what MVB++ builds 14_dab with
    auto mas = OpenMagneticsTesting::mas_loader(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'))
                                                + "/../MAS/examples/14_dab_xfmr_pm8770_n97.json");
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto layout = coil.get_connection_layout();
    size_t ramps = 0, straight = 0, stationWouldRamp = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
        const auto stub = pitchTrueStubOf(coil, route);
        INFO(route.winding << " parallel " << route.parallel
                           << (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? " entrance" : " exit")
                           << (stub ? " station h " + std::to_string(stub->station * 1e3) + " mm, at the slot " +
                                          std::to_string(stub->pitchTrue * 1e3) + " mm" : std::string(" (no stub)")));
        const double R = route.plannedBendRadius;
        if (!stub) {
            CHECK_FALSE(route.rampLength);
            continue;
        }
        // The two rules disagree on a lead whose stub is shorter than 2R at one of its station and its
        // slot and not at the other.
        if ((stub->station < 2 * R) != (stub->pitchTrue < 2 * R)) ++stationWouldRamp;
        if (route.rampLength) {
            CHECK(stub->pitchTrue < 2 * R);
            const double rho = (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.waypoints.back() : route.waypoints.front())[0];
            const double Rg = 1.0 / std::sqrt(1.0 / (R * R) - 1.0 / (rho * rho));   // PM: a round column
            CHECK_THAT(*route.rampLength, Catch::Matchers::WithinRel(std::numbers::pi * std::sqrt(stub->pitchTrue * Rg / 2), 1e-12));
            ++ramps;
        }
        else {
            CHECK(stub->pitchTrue >= 2 * R);
            ++straight;
        }
        // Spec S0/S1 (Alf, 2026-10-25): the Primary, the primary isolation group, packs on the x <= 0
        // half, so its p1 exit leaves at -7.635 mm, where its last turn stands 1.550 mm off its row
        // (b95d64cf: on the x >= 0 half, 3.77 mm). That is under the 2R = 1.655 mm its two bends need,
        // so it ramps: the stub ramps when it is shorter than 2R at its slot.
        if (route.winding == "Primary" && route.parallel == 1 && route.kind == ConnectionKind::TERMINAL_EXIT) {
            CHECK(stub->pitchTrue < 2 * R);
            CHECK(route.rampLength);
        }
        // And the other way round: the Primary's p2 exit stands 1.846 mm off its row at its station, 2R
        // clear, but 0.752 mm at its slot (-3.436 mm), so it ramps although the station rule would not.
        if (route.winding == "Primary" && route.parallel == 2 && route.kind == ConnectionKind::TERMINAL_EXIT) {
            CHECK(stub->station >= 2 * R);
            CHECK(stub->pitchTrue < 2 * R);
            CHECK(route.rampLength);
        }
    }
    UNSCOPED_INFO(ramps << " ramps, " << straight << " straight stubs, " << stationWouldRamp << " on which the station rule decides otherwise");
    REQUIRE(stationWouldRamp > 0);   // the premise: the two rules disagree on this design
    REQUIRE(ramps > 0);
    settings.reset();
}

// ABT #1336 (Alf, 2026-10-19: "MKF decides, MVB++ draws"): a consumer bends EVERY terminal lead at the
// radius MKF planned for it. A concentric terminal route carries it (plannedBendRadius); the toroidal
// emitter records no routes (its frame is polar), so MVB++ had nothing to read and refused every toroid.
// Coil::terminal_lead_bend_radius is the one rule both answer from.
TEST_CASE("Every terminal lead has MKF's planned bend radius, routed or not (ABT #1336)",
          "[constructive-model][coil][pins][abt1336]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const std::string here = std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/'));
    SECTION("a concentric coil: every terminal route carries exactly the rule's radius") {
        settings.reset();
        clear_databases();
        settings.set_coil_use_real_winding_geometry(true);
        auto mas = OpenMagneticsTesting::mas_loader(here + "/../MAS/examples/14_dab_xfmr_pm8770_n97.json");
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto& coil = magnetic.get_mutable_coil();
        REQUIRE(coil.is_real_winding_blocking_applied());
        size_t terminals = 0;
        for (const auto& route : coil.get_connection_layout().routes) {
            if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
            INFO(route.winding << " parallel " << route.parallel);
            const auto end = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? End::START : End::FINISH;
            CHECK(route.plannedBendRadius > 0);
            CHECK(route.plannedBendRadius == coil.terminal_lead_bend_radius(route.winding, end, route.parallel));
            ++terminals;
        }
        REQUIRE(terminals == 16);
        settings.reset();
    }
    SECTION("a toroid, whose terminal routes carry no planned bend: the rule still plans each lead's bend") {
        settings.reset();
        clear_databases();
        settings.set_coil_use_real_winding_geometry(true);
        auto mas = OpenMagneticsTesting::mas_loader(here + "/testData/toroidal_inductor_round_wire_multilayer.json");
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto& coil = magnetic.get_mutable_coil();
        // The premise (Alf, 2026-10-30, since toroid routes exist): the toroid's terminal routes carry no
        // plannedBendRadius, so the bend comes from the rule alone.
        size_t terminalRoutes = 0;
        for (const auto& route : coil.get_connection_layout().routes) {
            if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) continue;
            CHECK(route.toroidal);
            CHECK(route.plannedBendRadius == 0);
            ++terminalRoutes;
        }
        CHECK(terminalRoutes == 2);
        const auto wire = coil.resolve_wire(0);
        const double coatedRadius = std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height()) / 2;
        for (auto end : {End::START, End::FINISH}) {
            const double radius = coil.terminal_lead_bend_radius("Primary", end, 0);
            CHECK(radius >= coatedRadius);
            CHECK(radius == OpenMagnetics::Coil::lead_bend_radius(wire, std::nullopt, coatedRadius, coil.is_real_winding_blocking_applied()));
        }
        CHECK_THROWS(coil.terminal_lead_bend_radius("No such winding", End::START, 0));
        settings.reset();
    }
}
