// ABT #1705: MKF gives every terminal lead its own lane on the connection face.
//
// Fixture: 14_dab_xfmr_pm8770_n97 (PM 87/70, N97), Primary and Secondary each 12 turns of 4 parallels of
// Round 0.80 Grade 1 (coated OD 0.855 mm), as the C2 corpus processed it. Under real winding each
// winding is laid as an inner layer (11 turns per parallel) and an outer layer (the 12th turn,
// descending to the bottom flange), so a winding's ENTRANCE leads (inner layer) and EXIT leads (outer
// layer) leave by the same flange, on rows 0.855 mm apart: in the window section the two bundles
// merely touch, so the old lane rule let them share the face, and in 3D -- where every lead turns off
// its helix and onto its radial run through bends of its planned radius -- the inner-layer entrance
// leads ran out through the outer layer exactly where the exit leads attach (MVB++: 7 violations, the
// worst 815 um inside the coated envelope).
//
// Every expected value is derived here from the routes' own data (coated diameters from the wires,
// waypoints, plannedBendRadius, rampLength), independently of terminal_exit_slots.
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <map>
#include <numbers>
#include <optional>
#include <filesystem>
#include <fstream>
#include <limits>
#include <source_location>
#include "constructive_models/Wire.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/ResultsAnalyzer.h"
#include "support/Painter.h"
#include "TestingUtils.h"
#include "json.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <algorithm>
#include <numbers>
#include <iterator>
#include <sstream>

using json = nlohmann::json;
using OpenMagnetics::Settings;

namespace {

OpenMagnetics::Coil wind_dab_fixture() {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append("abt1705_dab_pm8770_4_parallels.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    // As MVB++ enriches it for the C2 corpus: real winding, corners planned at 1.05 coated radii
    // (MVB++ kRoundCornerBendFactor), no declared minimum bend, leads ending at the window border.
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    auto wound = enriched.get_coil();
    settings.reset();
    REQUIRE(wound.is_real_winding_blocking_applied());
    return wound;
}

}  // namespace

namespace {

using OpenMagnetics::ConnectionKind;
using OpenMagnetics::ConnectionRoute;
using OpenMagnetics::TerminalTurnPath;
using OpenMagnetics::LeadFinding;
using OpenMagnetics::LeadFindingKind;
using OpenMagnetics::TerminalLanePlan;
using WindingSides = std::vector<std::pair<std::string, IsolationSide>>;

bool is_terminal(const ConnectionRoute& route) {
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT;
}

std::string label(const ConnectionRoute& route) {
    return route.winding + " p" + std::to_string(route.parallel) +
           (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? " entrance" : " exit");
}

// The lane plus, for a ramped lead, the ramp and its radial bend's level run on the side its turn lies
// (ConnectionRoute::rampLength): an entrance's turn leaves towards -x, an exit's arrives from +x.
std::pair<double, double> footprint(const ConnectionRoute& route) {
    const double x = route.exitSlot.value();
    const double span = route.rampLength ? *route.rampLength + route.plannedBendRadius : 0.0;
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE ? std::pair{x - span, x} : std::pair{x, x + span};
}

double gap_between(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return std::max(a.first - b.second, b.first - a.second);
}

double segment_distance_2d(const std::vector<double>& p0, const std::vector<double>& p1, const std::vector<double>& q0,
                           const std::vector<double>& q1) {
    // Two segments in the plane: the minimum is at an endpoint of one against the other, or 0 when they cross.
    auto point_segment = [](const std::vector<double>& p, const std::vector<double>& a, const std::vector<double>& b) {
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        const double length2 = dx * dx + dy * dy;
        double t = length2 > 0 ? ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2 : 0.0;
        t = std::min(1.0, std::max(0.0, t));
        return std::hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
    };
    auto cross = [](const std::vector<double>& o, const std::vector<double>& a, const std::vector<double>& b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    };
    const double d1 = cross(q0, q1, p0);
    const double d2 = cross(q0, q1, p1);
    const double d3 = cross(p0, p1, q0);
    const double d4 = cross(p0, p1, q1);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
        return 0.0;
    }
    return std::min({point_segment(p0, q0, q1), point_segment(p1, q0, q1), point_segment(q0, p0, p1), point_segment(q1, p0, p1)});
}

double section_distance(const ConnectionRoute& a, const ConnectionRoute& b) {
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < a.waypoints.size(); ++i) {
        for (size_t j = 0; j + 1 < b.waypoints.size(); ++j) {
            best = std::min(best, segment_distance_2d(a.waypoints[i], a.waypoints[i + 1], b.waypoints[j], b.waypoints[j + 1]));
        }
    }
    return best;
}

// ABT #1336: a round column's turn paths, from each terminal route's attach waypoint radius
// (Coil::terminal_turn_path's round branch, for hand-built routes that have no coil).
// Hand-built routes are unsleeved: each turn and stub is the route's own diameter.
std::vector<std::optional<TerminalTurnPath>> round_paths(const std::vector<ConnectionRoute>& routes, const std::vector<double>& diameters) {
    REQUIRE(diameters.size() == routes.size());
    std::vector<std::optional<TerminalTurnPath>> paths;
    for (size_t index = 0; index < routes.size(); ++index) {
        const auto& route = routes[index];
        const auto& attach = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.waypoints.back() : route.waypoints.front();
        TerminalTurnPath path{ColumnShape::ROUND, 2 * std::numbers::pi * attach.at(0), 0.0, attach.at(0), attach.at(0)};
        path.turnWireDiameter = diameters[index];
        path.stubDiameter = diameters[index];
        paths.push_back(path);
    }
    return paths;
}

}  // namespace

TEST_CASE("14_dab: every terminal lead has its own lane, clear of every other bundle, with room for its bends (ABT #1705)",
          "[constructive-model][coil][real-winding][abt-1705]") {
    auto coil = wind_dab_fixture();
    const auto layout = coil.get_connection_layout();
    const auto wires = coil.get_wires();
    std::vector<ConnectionRoute> leads;
    for (const auto& route : layout.routes) {
        if (is_terminal(route)) {
            leads.push_back(route);
        }
    }
    // 2 windings x 4 parallels x (entrance, exit), every one with a lane (a round column: no tangential departures).
    REQUIRE(leads.size() == 16);
    auto diameter_of = [&](const ConnectionRoute& route) {
        auto wire = wires[coil.get_winding_index_by_name(route.winding)];
        return std::max({wire.get_maximum_outer_width(), wire.get_maximum_outer_height(), route.sleeveOuterDiameter.value_or(0.0)});
    };
    for (const auto& lead : leads) {
        INFO(label(lead));
        REQUIRE(lead.exitSlot.has_value());
        REQUIRE_FALSE(lead.tangentDeparture.has_value());
        REQUIRE(lead.plannedBendRadius >= diameter_of(lead) / 2);
    }

    const auto coilTurns = coil.get_turns_description().value();
    auto path_of = [&](const ConnectionRoute& route) {
        const std::string& name = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.toTurn : route.fromTurn;
        const auto turn = std::find_if(coilTurns.begin(), coilTurns.end(), [&](const auto& candidate) { return candidate.get_name() == name; });
        REQUIRE(turn != coilTurns.end());
        return coil.terminal_turn_path(route, *turn, layout.routes);
    };
    size_t crossBundlePairs = 0;
    size_t primaryEntranceVsExit = 0;
    double minimumCrossBundleSurplus = std::numeric_limits<double>::max();
    double minimumCrossBundleSlack = std::numeric_limits<double>::max();   // gap beyond the required one
    for (size_t i = 0; i < leads.size(); ++i) {
        for (size_t j = i + 1; j < leads.size(); ++j) {
            const auto& a = leads[i];
            const auto& b = leads[j];
            if (a.side != b.side) {
                continue;   // another isolation side's face
            }
            INFO(label(a) << " at x " << *a.exitSlot * 1e3 << " mm vs " << label(b) << " at x " << *b.exitSlot * 1e3 << " mm");
            const double radii = (diameter_of(a) + diameter_of(b)) / 2;
            const double gap = gap_between(footprint(a), footprint(b));
            const double section = section_distance(a, b);
            // No two leads within their summed coated radii, siblings included.
            CHECK((gap >= radii - 1e-9 || section >= radii - 1e-9));
            const bool sameBundle = a.winding == b.winding && a.kind == b.kind;
            if (sameBundle) {
                continue;
            }
            // Two bundles whose leads can meet in 3D stand apart along the face: disjoint, with a clearance
            // gap, never an exact touch. Alf, 2026-10-28: "can meet" is the lane planner's own 3D judgement
            // (spec P5: real sections, bend arcs and the attach turns' spans), not the section-distance
            // proxy this test used before (routes within the coated radii plus both bends in the section),
            // which demanded face room between Secondary p1's entrance ramp and p2's exit, about 1.1 mm
            // apart at their coated surfaces in 3D.
            if (OpenMagnetics::Coil::terminal_leads_meet_as_planned(a, path_of(a), b, path_of(b))) {
                ++crossBundlePairs;
                if (a.winding == "Primary" && b.winding == "Primary") {
                    ++primaryEntranceVsExit;
                }
                // Alf 2026-10-25 (C3): no blanket bend-radius gap between bundles; S1 puts the entrance on -x.
                // Two bundles stand their coated radii apart along the face (14_dab: the coated OD, 0.855 mm).
                const double required = radii;
                CHECK(gap >= required - 1e-9);
                minimumCrossBundleSurplus = std::min(minimumCrossBundleSurplus, gap - radii);
                minimumCrossBundleSlack = std::min(minimumCrossBundleSlack, gap - required);
            }
        }
    }
    // The pairs this test is about exist: the Primary's inner-layer entrances against its outer-layer exits.
    REQUIRE(primaryEntranceVsExit == 16);
    REQUIRE(crossBundlePairs >= primaryEntranceVsExit);
    INFO("smallest cross-bundle clearance beyond the coated radii: " << minimumCrossBundleSurplus * 1e3 << " mm");
    // Alf 2026-10-25 (C3): no blanket bend-radius gap between bundles; the surplus beyond the coated radii is 0.
    CHECK(minimumCrossBundleSurplus > -1e-9);
    // The tightest bundle boundary stands EXACTLY the required gap apart: the anchor is solved, not
    // stepped along a pitch grid that would waste up to a pitch at every boundary (Alf, ABT #1705).
    INFO("smallest cross-bundle slack beyond the required gap: " << minimumCrossBundleSlack * 1e3 << " mm");
    CHECK(std::abs(minimumCrossBundleSlack) < 1e-9);

    // Each lead leaves its helix with room for its bends: a straight stub from its turn to its row
    // holds two bends of plannedBendRadius (one in the turn surface, one onto the radial run), so it is
    // at least 2R long, or it is published as a ramp on the turn surface (ABT #1336) whose span the
    // lanes above already keep clear. The stub is measured PITCH-TRUE (ABT #1336, Alf, 2026-10-19):
    // from where its own turn really stands at its lane (Coil::terminal_turn_height_at), not from the
    // turn's 2-D station at the connection plane.
    const auto turns = coil.get_turns_description().value();
    auto attach_advance = [&](const ConnectionRoute& lead, const MAS::Turn*& attachTurn) {
        const bool entrance = lead.kind == ConnectionKind::TERMINAL_ENTRANCE;
        const std::string& turnName = entrance ? lead.toTurn : lead.fromTurn;
        size_t at = turns.size();
        for (size_t k = 0; k < turns.size(); ++k) {
            if (turns[k].get_name() == turnName) at = k;
        }
        REQUIRE(at < turns.size());
        attachTurn = &turns[at];
        auto same = [&](size_t k) {
            return turns[k].get_winding() == turns[at].get_winding() && turns[k].get_parallel() == turns[at].get_parallel();
        };
        if (entrance) {
            for (size_t k = at + 1; k < turns.size(); ++k) {
                if (same(k)) return turns[k].get_coordinates()[1] - turns[at].get_coordinates()[1];
            }
        }
        else {
            for (size_t k = at; k-- > 0;) {
                if (same(k)) return turns[at].get_coordinates()[1] - turns[k].get_coordinates()[1];
            }
        }
        return 0.0;
    };
    size_t stubs = 0;
    for (const auto& lead : leads) {
        // The route runs from the turn (entrance: last waypoint; exit: first) to the border.
        const auto& turnEnd = lead.kind == ConnectionKind::TERMINAL_ENTRANCE ? lead.waypoints.back() : lead.waypoints.front();
        const auto& next = lead.kind == ConnectionKind::TERMINAL_ENTRANCE ? lead.waypoints[lead.waypoints.size() - 2] : lead.waypoints[1];
        const bool axialStub = std::abs(turnEnd[0] - next[0]) < 1e-12 && std::abs(turnEnd[1] - next[1]) > 0;
        if (!axialStub) {
            continue;
        }
        ++stubs;
        const MAS::Turn* attachTurn = nullptr;
        const double advance = attach_advance(lead, attachTurn);
        const auto path = coil.terminal_turn_path(lead, *attachTurn, layout.routes);
        const double turnY = OpenMagnetics::Coil::terminal_turn_height_at(lead, advance, path, lead.exitSlot.value());
        const double height = std::abs(turnY - next[1]);
        INFO(label(lead) << ": stub " << height * 1e3 << " mm, bend radius " << lead.plannedBendRadius * 1e3 << " mm");
        CHECK((height >= 2 * lead.plannedBendRadius - 1e-12 || lead.rampLength.has_value()));
    }
    REQUIRE(stubs > 0);
}

TEST_CASE("Pins decide the side: with no pins the S0/S1 default stands, one face, halves by isolation side (spec section 12, ABT #1336)",
          "[constructive-model][coil][real-winding][abt-1705][abt1336]") {
    auto coil = wind_dab_fixture();
    // (iii) no pins: the premise of this case.
    const auto bobbin = coil.resolve_bobbin();
    REQUIRE(bobbin.get_processed_description());
    REQUIRE((!bobbin.get_processed_description()->get_pins() || bobbin.get_processed_description()->get_pins()->empty()));
    const auto layout = coil.get_connection_layout();
    std::set<int> halves;
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        INFO(label(route));
        CHECK(route.side == 0);   // S0: ONE FACE
        if (!is_terminal(route)) {
            continue;
        }
        const auto& winding = coil.get_functional_description()[coil.get_winding_index_by_name(route.winding)];
        REQUIRE(route.laneHalf.has_value());
        CHECK(route.laneHalf.value() == (winding.get_isolation_side() == MAS::IsolationSide::PRIMARY ? -1 : 1));   // S0 halves
        halves.insert(route.laneHalf.value());
        ++terminals;
    }
    REQUIRE(terminals == 16);
    REQUIRE(halves.size() == 2);   // both isolation groups present
}

TEST_CASE("check_terminal_lanes refuses lanes that put two bundles' bent leads together (ABT #1705)",
          "[constructive-model][coil][real-winding][abt-1705]") {
    auto coil = wind_dab_fixture();
    const auto routes = coil.get_connection_layout().routes;
    const auto wires = coil.get_wires();
    std::vector<double> diameters;
    std::vector<std::optional<double>> slots;
    std::optional<size_t> entrance;
    std::optional<size_t> exit;
    for (size_t index = 0; index < routes.size(); ++index) {
        auto wire = wires[coil.get_winding_index_by_name(routes[index].winding)];
        diameters.push_back(std::max({wire.get_maximum_outer_width(), wire.get_maximum_outer_height()}));
        slots.push_back(routes[index].exitSlot);
        if (routes[index].winding == "Primary" && routes[index].parallel == 0 && routes[index].kind == ConnectionKind::TERMINAL_ENTRANCE) {
            entrance = index;
        }
        if (routes[index].winding == "Primary" && routes[index].parallel == 3 && routes[index].kind == ConnectionKind::TERMINAL_EXIT) {
            exit = index;
        }
    }
    REQUIRE(entrance.has_value());
    REQUIRE(exit.has_value());
    // Each terminal route's attach turn path (entrance: its toTurn, exit: its fromTurn); none for the others.
    const auto turns = coil.get_turns_description().value();
    std::vector<std::optional<TerminalTurnPath>> turnPaths;
    for (const auto& route : routes) {
        if (!is_terminal(route)) {
            turnPaths.push_back(std::nullopt);
            continue;
        }
        const std::string& turnName = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.toTurn : route.fromTurn;
        const MAS::Turn* attachTurn = nullptr;
        for (const auto& turn : turns) {
            if (turn.get_name() == turnName) attachTurn = &turn;
        }
        REQUIRE(attachTurn != nullptr);
        turnPaths.push_back(coil.terminal_turn_path(route, *attachTurn, routes));
    }
    // MKF's own lanes pass: no lead carries a finding.
    std::vector<std::vector<LeadFinding>> findings;
    REQUIRE_NOTHROW(findings = OpenMagnetics::Coil::check_terminal_lanes(routes, diameters, slots, true, turnPaths));
    REQUIRE(findings.size() == routes.size());
    for (const auto& routeFindings : findings) {
        CHECK(routeFindings.empty());
    }
    // The lanes the old rule handed out put the Primary's parallel 3 exit on its parallel 0 entrance's lane:
    // their rows touch in the window section (0.855 mm), so the old test saw no clash.
    // Spec section 17 (Alf 2026-10-26): under real winding the collision is reported on both leads, not
    // thrown (was CHECK_THROWS_WITH the same message).
    slots[exit.value()] = slots[entrance.value()];
    REQUIRE_NOTHROW(findings = OpenMagnetics::Coil::check_terminal_lanes(routes, diameters, slots, true, turnPaths));
    for (size_t lead : {exit.value(), entrance.value()}) {
        REQUIRE_FALSE(findings[lead].empty());
        CHECK(findings[lead].front().kind == LeadFindingKind::LANES_COLLIDE);
        CHECK_THAT(findings[lead].front().witness, Catch::Matchers::ContainsSubstring("Terminal lanes collide") &&
                                                      Catch::Matchers::ContainsSubstring("Primary' parallel 3 exit") &&
                                                      Catch::Matchers::ContainsSubstring("different bundles"));
    }
    // Two leads on one point collide even without the bend rule.
    auto coincident = routes;
    coincident[exit.value()].waypoints = coincident[entrance.value()].waypoints;
    CHECK_THROWS_WITH(OpenMagnetics::Coil::check_terminal_lanes(coincident, diameters, slots, false, {}),
                      Catch::Matchers::ContainsSubstring("their coated radii need"));
}

TEST_CASE("terminal_exit_slots: two bundles that touch in the section take disjoint lanes only under real winding (ABT #1705)",
          "[constructive-model][coil][abt-1705]") {
    // 14_dab's Primary in miniature (mm -> m): the inner-layer entrance row and the outer-layer exit row,
    // one coated OD apart in the section; two parallels each, side by side on their rows.
    const double od = 0.855e-3;
    const double bend = 1.05 * od / 2;
    auto route = [&](int64_t parallel, ConnectionKind kind, std::vector<std::vector<double>> waypoints) {
        ConnectionRoute r;
        r.winding = "Primary";
        r.parallel = parallel;
        r.kind = kind;
        r.side = 0;
        r.waypoints = waypoints;
        r.plannedBendRadius = bend;
        (kind == ConnectionKind::TERMINAL_ENTRANCE ? r.toTurn : r.fromTurn) = "t";
        // Spec section 21: 14_dab's turns are wound clockwise (a window right of the column), so each
        // attach turn leaves the crossing towards -x (Coil::terminal_wrap_direction).
        r.wrapDirection = -1;
        return r;
    };
    std::vector<ConnectionRoute> routes = {
        route(0, ConnectionKind::TERMINAL_ENTRANCE, {{34.5025e-3, -21.4223e-3}, {18.8475e-3, -21.4223e-3}}),
        route(1, ConnectionKind::TERMINAL_ENTRANCE, {{34.5025e-3, -21.4223e-3}, {18.8475e-3, -21.4223e-3}, {18.8475e-3, -19.6957e-3}}),
        route(0, ConnectionKind::TERMINAL_EXIT, {{19.7025e-3, -18.7212e-3}, {19.7025e-3, -20.5673e-3}, {34.5025e-3, -20.5673e-3}}),
        route(1, ConnectionKind::TERMINAL_EXIT, {{19.7025e-3, -17.8086e-3}, {19.7025e-3, -20.5673e-3}, {34.5025e-3, -20.5673e-3}}),
    };
    const std::vector<double> attach = {-21.4223e-3, -19.6957e-3, -18.7212e-3, -17.8086e-3};
    const std::vector<double> diameters(routes.size(), od);
    // The stubs differ in every bundle here, so the attach advance never decides an order.
    const std::vector<double> advance(routes.size(), 3.6618e-3);
    const WindingSides primaryOnly = {{"Primary", IsolationSide::PRIMARY}};

    SECTION("ideal winding: the rows touch, so the old rule lets both bundles start at the plane") {
        auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, false, advance, {}, primaryOnly);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == od);
        CHECK(slots[3].value() == 0.0);   // exit order by span descending: parallel 1 first
        CHECK(slots[2].value() == od);
    }
    SECTION("real winding: the exit bundle starts exactly d + R past the entrance bundle") {
        auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, true, advance, round_paths(routes, diameters), primaryOnly);
        // Spec section 25 (Alf, approved 2026-10-28): the entrance order comes from the wrap directions.
        // Parallel 1's stub climbs past parallel 0's first turn, which leaves its lane towards -x, so
        // parallel 1 stands at the crossing and parallel 0 one OD out (was p0 0, p1 -od by stub span).
        CHECK(slots[1].value() == 0.0);
        // Alf 2026-10-25 (C3): no blanket bend-radius gap between bundles; S1 puts the entrance on -x.
        CHECK(slots[0].value() == -od);
        // The exits start one coated OD past the crossing (0.855 mm), then one OD apart (1.710 mm).
        const double firstExit = od;
        CHECK(std::abs(slots[3].value() - firstExit) < 1e-12);
        CHECK(std::abs(slots[2].value() - (firstExit + od)) < 1e-12);
    }
    SECTION("real winding without a planned bend throws") {
        routes[2].plannedBendRadius = 0;
        CHECK_THROWS_WITH(OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, true, advance, round_paths(routes, diameters), primaryOnly),
                          Catch::Matchers::ContainsSubstring("planned bend radius"));
    }
}

// ============================================================================================
// ABT #1383: the terminal-lead exit of CLOSED cores (drum + ring SEATED on a base plate; molded
// body with a lead frame), planned by Coil::plan_closed_core_lead_exits and carried on
// ConnectionLayout (routes' exitRun, leadExitSolids) for MVB++ (#1381) to draw.
//
// Fixtures (generic models; each carries its provenance in magnetic.manufacturerInfo.datasheetInfo,
// saying which values are published and which are inferred or stand-ins):
// - lead_exit_drumring_7332.json: a 7.3 x 7.3 x 3.2 mm drum + ring part, the ring SEATED on a square
//   base plate M = 7.3 mm (the footprint), N = 0 (no chamfer value published), B = the published
//   3.2 mm part height with E = B - D - F and L = B - F, pads 1.5 x 1.5 mm; NiZn ferrite (the
//   published material class; the grade is a stand-in); wire inferred from the published DCR.
// - lead_exit_molded_4030.json: a 4.0 x 4.0 x 3.1 mm molded part, pads 3.25 x 0.82 mm; wire inferred
//   from the published DCR; the powder material is a stand-in.
// - lead_exit_molded_overfull.json: test geometry, not a part: 7 turns of 0.355 mm wire in the same
//   cavity; the copper fits but the lead lane passes the body side, which keeps draw-and-report tested.
// Both bobbins 'Dummy': a drum is its own former and a molded coil has none.
// (Moved here from the lead-exit branch's own test file: terminal-lead tests live together.)
namespace lead_exit_suite {

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::StartsWith;

namespace {

json load_fixture(const std::string& name) {
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), name));
    REQUIRE(file.good());
    return json::parse(file);
}

double dimension(const json& fixture, const std::string& letter) {
    return fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"][letter]["nominal"].get<double>();
}

OpenMagnetics::Magnetic autocomplete(json magneticJson) {
    OpenMagnetics::Magnetic magnetic(magneticJson);
    return magnetic_autocomplete(magnetic);
}

double polyline_length(const std::vector<std::vector<double>>& points, size_t from = 0) {
    double length = 0;
    for (size_t k = from; k + 1 < points.size(); ++k) {
        double sum = 0;
        for (size_t axis = 0; axis < points[k].size(); ++axis) {
            sum += std::pow(points[k + 1][axis] - points[k][axis], 2);
        }
        length += std::sqrt(sum);
    }
    return length;
}

bool is_terminal(const ConnectionRoute& route) {
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT;
}

size_t count_kind(const ConnectionLayout& layout, LeadExitSolid::Kind kind) {
    size_t count = 0;
    for (const auto& solid : layout.leadExitSolids) {
        count += solid.kind == kind;
    }
    return count;
}

// Every terminal route's routedLength is its in-window run plus its exit run (arc at arc length).
void check_routed_lengths(const ConnectionLayout& layout) {
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        REQUIRE(route.exitRunArc);
        const auto& arc = route.exitRunArc.value();
        const double expected = polyline_length(route.waypoints) + arc.radius * std::abs(arc.toAzimuth - arc.fromAzimuth) +
                                polyline_length(route.exitRun, 1);
        CHECK_THAT(route.routedLength, WithinAbs(expected, 1e-8));
    }
}

}  // namespace

TEST_CASE("A drum with its ring SEATED on a square base plate lays its leads out through ring notches to its electrodes",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::RING_NOTCH);

    const auto layout = coil.get_connection_layout();
    const auto spaces = coil.get_connection_reserved_spaces();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    const double plateTop = -dimension(fixture, "E") / 2;   // the window floor: groove E, centred
    const double ringOuterRadius = dimension(fixture, "J") / 2;

    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        // An edge-routed lead is on the plate-side (bottom) row: its in-window run's outer end is
        // below the window centre. (An outermost-layer end runs radially at its own turn's height.)
        const auto& border = entrance ? route.waypoints.front() : route.waypoints.back();
        for (const auto& space : spaces) {
            if (space.isTerminal && space.kind == route.kind && space.edgeDepth > 0) {
                CHECK(border[1] < 0);
            }
        }
        REQUIRE(route.exitRun.size() >= 3);
        // Every lead passes its notch on the plate-side rows: the exit run's point just outside the
        // ring (the start of its corner run) is below the window centre and on or above the plate.
        std::optional<std::vector<double>> found;
        for (const auto& point : route.exitRun) {
            if (std::abs(std::hypot(point[0], point[2]) - (ringOuterRadius + diameter / 2)) < 1e-9) {
                found = point;
            }
        }
        REQUIRE(found);
        const auto& throughNotch = found.value();
        CHECK(throughNotch[1] < 0);
        CHECK(throughNotch[1] - diameter / 2 >= plateTop - 1e-9);
        REQUIRE(route.exitRunEnd);
        CHECK(route.exitRunEnd.value() == ConnectionRoute::ExitRunEnd::ELECTRODE);
        // The run ends over its own side's terminal: entrance -X, exit +X, both on the -Z half.
        const auto& joint = route.exitRun.back();
        CHECK((entrance ? joint[0] < -ringOuterRadius : joint[0] > ringOuterRadius));
        CHECK(joint[2] < 0);
        // ... lying on the plate top.
        CHECK_THAT(joint[1], WithinAbs(plateTop + diameter / 2, 1e-9));
    }
    CHECK(terminals == 2);
    check_routed_lengths(layout);

    // Two notches, on the diagonals to the two connection-plane corners.
    REQUIRE(count_kind(layout, LeadExitSolid::Kind::CORE_CUT) == 2);
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind != LeadExitSolid::Kind::CORE_CUT) {
            continue;
        }
        REQUIRE(solid.body.xz.size() == 4);
        double cx = 0, cz = 0;
        for (const auto& point : solid.body.xz) {
            cx += point[0] / 4;
            cz += point[1] / 4;
        }
        CHECK(cz < 0);
        CHECK(solid.entrance == (cx < 0));
        CHECK_THAT(std::abs(cx), WithinAbs(std::abs(cz), 1e-9));   // on the 45 degree diagonal
        // Width = parallels (1) x the largest coated OD: the inner edge xz[0] -> xz[3].
        const double width = std::hypot(solid.body.xz[3][0] - solid.body.xz[0][0], solid.body.xz[3][1] - solid.body.xz[0][1]);
        CHECK_THAT(width, WithinAbs(1 * diameter, 1e-9));
        // Through the ring wall: the outer edge at the ring's outer radius J/2.
        const double outerCentre = std::hypot((solid.body.xz[1][0] + solid.body.xz[2][0]) / 2, (solid.body.xz[1][1] + solid.body.xz[2][1]) / 2);
        CHECK_THAT(outerCentre, WithinAbs(ringOuterRadius, 1e-9));
        // Height = the lead's height above the plate: its edge row's edgeDepth, or one OD for an
        // outermost-layer end, which descends to lie on the plate.
        CHECK_THAT(solid.body.y0, WithinAbs(plateTop, 1e-12));
        std::optional<double> expectedHeight;
        for (const auto& space : spaces) {
            if (!space.isTerminal || space.kind != (solid.entrance ? ConnectionKind::TERMINAL_ENTRANCE : ConnectionKind::TERMINAL_EXIT)) {
                continue;
            }
            if (space.edgeDepth > 0) {
                expectedHeight = std::max(expectedHeight.value_or(0.0), space.edgeDepth);
            }
        }
        CAPTURE(solid.name, solid.body.y1 - solid.body.y0);
        CHECK_THAT(solid.body.y1 - solid.body.y0, WithinAbs(expectedHeight.value_or(diameter), 1e-9));
        CHECK(solid.material.empty());
    }
    // Each end's terminal: a side leg and a mounting pad (the board face), of the sourced stock.
    CHECK(count_kind(layout, LeadExitSolid::Kind::ELECTRODE) == 4);
    size_t boardFaces = 0;
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind == LeadExitSolid::Kind::ELECTRODE) {
            CHECK(solid.material == "phosphor bronze");
            boardFaces += solid.boardFace.has_value();
        }
    }
    CHECK(boardFaces == 2);
    settings.reset();
}

TEST_CASE("A bifilar drum-ring end leaves through one notch, its parallels side by side, every row on the plate",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    // Two leads per end side by side on the 1.5 mm terminal. With the fixture's wire both the ideal
    // and the real winding leave room for the two lanes inside the ring bore (the refusal when they
    // do not is covered in "A closed core's lead exit refuses what it cannot build").
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    const auto layout = coil.get_connection_layout();
    const auto spaces = coil.get_connection_reserved_spaces();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    // Every edge-routed terminal row is on the plate-side (bottom) edge.
    size_t edgeRows = 0;
    for (const auto& space : spaces) {
        if (space.isTerminal && space.edgeDepth > 0 &&
            (space.kind == ConnectionKind::TERMINAL_ENTRANCE || space.kind == ConnectionKind::TERMINAL_EXIT)) {
            ++edgeRows;
            CAPTURE(space.coordinates[1]);
            CHECK(space.coordinates[1] < 0);
        }
    }
    CHECK(edgeRows > 0);
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        terminals += is_terminal(route) && !route.exitRun.empty();
    }
    CHECK(terminals == 4);
    // One notch per END, two parallels wide.
    REQUIRE(count_kind(layout, LeadExitSolid::Kind::CORE_CUT) == 2);
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind != LeadExitSolid::Kind::CORE_CUT) {
            continue;
        }
        CHECK_FALSE(solid.parallel.has_value());
        const double width = std::hypot(solid.body.xz[3][0] - solid.body.xz[0][0], solid.body.xz[3][1] - solid.body.xz[0][1]);
        CHECK_THAT(width, WithinAbs(2 * diameter, 1e-9));
    }
    check_routed_lengths(layout);
    settings.reset();
}

TEST_CASE("A closed core's terminal leads are laid on the plate-side edge even where the top edge is nearer",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    // Four parallels: the first layer holds the bundle's first turns in slots 1-4 of 6, so the
    // fourth parallel's entrance turn sits ABOVE the window centre and its nearest edge is the top.
    // A closed core's leads leave over its base, so that entrance must still be laid on the bottom.
    // The ring is widened (K 8.0, J = M 9.6 mm) and the pad made 2.0 mm wide so the four lanes and
    // the four leads side by side fit: test geometry, not a catalogue part.
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    auto& winding = fixture["magnetic"]["coil"]["functionalDescription"][0];
    winding["numberTurns"] = 3;
    winding["numberParallels"] = 4;
    for (auto& connection : winding["connections"]) {
        connection["padWidth"] = 0.002;
    }
    auto& dimensions = fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"];
    dimensions["K"]["nominal"] = 0.008;
    dimensions["J"]["nominal"] = 0.0096;
    dimensions["M"]["nominal"] = 0.0096;
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    // The premise: some entrance turn is above the window centre.
    bool entranceAboveCentre = false;
    const auto layout = coil.get_connection_layout();
    const auto turns = coil.get_turns_description().value();
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE) {
            continue;
        }
        for (const auto& turn : turns) {
            if (turn.get_name() == route.toTurn && turn.get_coordinates()[1] > 0) {
                entranceAboveCentre = true;
            }
        }
    }
    REQUIRE(entranceAboveCentre);
    size_t edgeRows = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.isTerminal && space.edgeDepth > 0 &&
            (space.kind == ConnectionKind::TERMINAL_ENTRANCE || space.kind == ConnectionKind::TERMINAL_EXIT)) {
            ++edgeRows;
            CAPTURE(space.coordinates[1]);
            CHECK(space.coordinates[1] < 0);
        }
    }
    CHECK(edgeRows >= 4);
    size_t exitRuns = 0;
    for (const auto& route : layout.routes) {
        exitRuns += is_terminal(route) && !route.exitRun.empty();
    }
    CHECK(exitRuns == 8);
    settings.reset();
}

TEST_CASE("A molded body's leads are welded to lead frames under the cavity floor",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_molded_4030.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::LEAD_FRAME);
    // The winding fits its cavity (ideal and real winding alike), so the lead exit is laid out.
    double outerRadius = 0;
    const auto layers = coil.get_layers_description().value();
    for (const auto& layer : layers) {
        outerRadius = std::max(outerRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    }
    CAPTURE(outerRadius);
    CHECK(outerRadius <= dimension(fixture, "E") / 2);
    const auto layout = coil.get_connection_layout();
    const auto spaces = coil.get_connection_reserved_spaces();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    const double frameThickness = 0.000127;   // src/data/dfm_rules.json leadExitStock
    const double bodyWidth = dimension(fixture, "A"), bodyDepth = dimension(fixture, "C");
    const double bodyHeight = dimension(fixture, "B"), cavityHeight = dimension(fixture, "D");

    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        // An edge-routed lead is on the plate-side (bottom) row (an outermost-layer end runs
        // radially at its own turn's height and descends at its tab).
        const auto& border = entrance ? route.waypoints.front() : route.waypoints.back();
        for (const auto& space : spaces) {
            if (space.isTerminal && space.kind == route.kind && space.edgeDepth > 0) {
                CHECK(border[1] < 0);
            }
        }
        REQUIRE(route.exitRunEnd);
        CHECK(route.exitRunEnd.value() == ConnectionRoute::ExitRunEnd::LEAD_FRAME_WELD);
        // Ends lying on the cavity floor, over its weld.
        CHECK_THAT(route.exitRun.back()[1], WithinAbs(-cavityHeight / 2 + diameter / 2, 1e-9));
    }
    CHECK(terminals == 2);
    check_routed_lengths(layout);

    REQUIRE(count_kind(layout, LeadExitSolid::Kind::LEAD_FRAME) == 2);
    REQUIRE(count_kind(layout, LeadExitSolid::Kind::WELD) == 2);
    double entranceInner = std::numeric_limits<double>::lowest();
    double exitInner = std::numeric_limits<double>::max();
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind == LeadExitSolid::Kind::LEAD_FRAME) {
            CHECK(solid.material == "phosphor bronze");
            // In the frame plane under the cavity floor, inside the base pressing and the footprint.
            CHECK_THAT(solid.body.y1, WithinAbs(-cavityHeight / 2, 1e-12));
            CHECK_THAT(solid.body.y1 - solid.body.y0, WithinAbs(frameThickness, 1e-12));
            CHECK(solid.body.y0 >= -bodyHeight / 2 - 1e-12);
            for (const auto& point : solid.body.xz) {
                CHECK(std::abs(point[0]) <= bodyWidth / 2 + 1e-12);
                CHECK(std::abs(point[1]) <= bodyDepth / 2 + 1e-12);
                if (solid.entrance) {
                    entranceInner = std::max(entranceInner, point[0]);
                }
                else {
                    exitInner = std::min(exitInner, point[0]);
                }
            }
        }
        if (solid.kind == LeadExitSolid::Kind::WELD) {
            // One coated OD long and wide, from the frame to the top of the wire lying on it.
            CHECK_THAT(solid.body.y0, WithinAbs(-cavityHeight / 2, 1e-12));
            CHECK_THAT(solid.body.y1 - solid.body.y0, WithinAbs(diameter, 1e-12));
            CHECK_THAT(solid.body.xz[1][0] - solid.body.xz[0][0], WithinAbs(diameter, 1e-12));
            CHECK_THAT(solid.body.xz[2][1] - solid.body.xz[1][1], WithinAbs(diameter, 1e-12));
        }
    }
    // The two frames stand at least one stock thickness apart.
    CAPTURE(entranceInner, exitInner);
    CHECK(exitInner - entranceInner >= frameThickness - 1e-12);
    settings.reset();
}

TEST_CASE("A molded body whose lead lane passes its side draws the lead and reports it",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_molded_overfull.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::LEAD_FRAME);
    // Real winding lays these 7 turns of 0.355 mm in TWO layers whose outer face (1.567 mm) is
    // inside the 1.75 mm cavity: the copper fits the body...
    double outerRadius = 0;
    const auto layers = coil.get_layers_description().value();
    REQUIRE(layers.size() == 2);
    for (const auto& layer : layers) {
        outerRadius = std::max(outerRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
        CHECK(std::abs(layer.get_coordinates()[1]) + layer.get_dimensions()[1] / 2 <= dimension(fixture, "D") / 2 + 1e-12);
    }
    CAPTURE(outerRadius);
    CHECK(outerRadius < dimension(fixture, "E") / 2);
    // ... but its lead lane does not: outer face + 0.384 mm ride-over + half the 0.384 mm pitch puts
    // the lane at 2.142 mm, 2.334 mm with its coated radius, past the body's side at 2.050 mm.
    // Drawn and reported, never refused: each lead is laid where the plan puts it, welded to its
    // lead frame, and the layout records one collision per lead.
    const auto layout = coil.get_connection_layout();
    size_t pastSide = 0;
    for (const auto& clash : layout.leadExitClashes) {
        if (clash.kind == LeadExitClash::Kind::LANE_PAST_BODY_SIDE) {
            ++pastSide;
            CHECK_THAT(clash.description, ContainsSubstring("past the body's side"));
        }
    }
    CHECK(pastSide == 2);
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (is_terminal(route)) {
            ++terminals;
            CHECK(route.exitRunEnd == ConnectionRoute::ExitRunEnd::LEAD_FRAME_WELD);
        }
    }
    CHECK(terminals == 2);
    // ... and analyze_results says so, once per clash.
    const auto findings = ResultsAnalyzer::analyze_lead_exit_clashes(magnetic);
    CHECK(findings.size() == layout.leadExitClashes.size());
    for (const auto& finding : findings) {
        CHECK(finding.code == ResultFindingCode::LEAD_EXIT_CLASHES);
        CHECK_THAT(finding.message, !ContainsSubstring("ABT"));
    }
    settings.reset();
}

// Alf, 2026-10-29: a closed core's lead-exit clashes are reported even when the turns' outlines cannot
// be judged (TURN_OUTLINES_NOT_EVALUABLE): the lead-exit plan records them without any turn outline.
// The molded body above, wound as it is, then its turns' stored shapes and its wire's outer diameter
// taken away, so that no turn has an outline.
TEST_CASE("A molded body's lead-exit clashes are reported when the turns' outlines are unknown",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_molded_overfull.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto coil = magnetic.get_coil();
    auto turns = coil.get_turns_description().value();
    for (auto& turn : turns) {
        turn.set_cross_sectional_shape(std::nullopt);
    }
    coil.set_turns_description(turns);
    auto wire = coil.resolve_wire(size_t(0));
    REQUIRE(wire.get_type() == WireType::ROUND);
    REQUIRE(wire.get_outer_diameter());
    wire.set_outer_diameter(std::nullopt);
    coil.get_mutable_functional_description()[0].set_wire(wire);
    magnetic.set_coil(coil);
    REQUIRE(magnetic.get_mutable_coil().turn_outlines_unknown_reason());
    const auto clashes = magnetic.get_mutable_coil().judge_lead_exit_clashes();
    REQUIRE_FALSE(clashes.empty());

    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    std::vector<ResultFinding> findings;
    REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(mas));
    size_t notEvaluable = 0;
    size_t reported = 0;
    for (const auto& finding : findings) {
        CAPTURE(finding.message);
        notEvaluable += finding.code == ResultFindingCode::TURN_OUTLINES_NOT_EVALUABLE;
        if (finding.code == ResultFindingCode::LEAD_EXIT_CLASHES) {
            ++reported;
            CHECK_THAT(finding.message, StartsWith("The ") && ContainsSubstring(" lead of turn ") &&
                                        !ContainsSubstring("Warning") && !ContainsSubstring("ABT"));
        }
    }
    CHECK(notEvaluable == 1);
    CHECK(reported == clashes.size());
    settings.reset();
}

TEST_CASE("The fixtures' inferred windings reproduce their published DC resistance",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    // The wire of each fixture is INFERRED from its published DC resistance (see its provenance);
    // this pins that inference: under real winding (turns as laid plus the routed terminal leads)
    // the MKF DC resistance at 25 C lies within the published typical-to-maximum spread of the
    // typical value. The inductance is deliberately NOT checked against the published one.
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    const std::string name = GENERATE(std::string("lead_exit_molded_4030.json"), std::string("lead_exit_drumring_7332.json"));
    CAPTURE(name);
    auto fixture = load_fixture(name);
    const auto& published = fixture["magnetic"]["manufacturerInfo"]["datasheetInfo"]["electrical"][0]["dcResistance"];
    const double typical = published["nominal"].get<double>();
    const double maximum = published["maximum"].get<double>();
    auto magnetic = autocomplete(fixture["magnetic"]);
    const double resistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_mutable_coil(), 25)[0];
    CAPTURE(resistance, typical, maximum);
    CHECK_THAT(resistance, WithinAbs(typical, maximum - typical));
    settings.reset();
}

TEST_CASE("A catalogue drumRing (no base plate M) has no lead exit, and says so",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["core"]["functionalDescription"]["shape"] = "DR 5.3x3.6 + RIS 7.3x3.5x5.7";
    // Autocomplete and the magnetics stay as they were for a catalogue part...
    auto magnetic = autocomplete(fixture["magnetic"]);
    CHECK_FALSE(magnetic.get_mutable_coil().get_core_lead_exit());
    // ... and laying its leads out is refused, naming why.
    REQUIRE_THROWS_WITH(magnetic.get_mutable_coil().get_connection_layout(), ContainsSubstring("has no lead exit"));
    REQUIRE_THROWS_WITH(magnetic.get_mutable_core().get_lead_exit_geometry(), ContainsSubstring("ring-around-both-flanges"));
}

// Alf, 2026-10-26: laying out a catalogue drumRing's leads stays refused (get_connection_layout throws),
// but analyze_results checks the core for its base plate M BEFORE any lead judgement lays the leads out,
// and reports the part's leads as not modelled instead of throwing.
TEST_CASE("A catalogue drumRing under real winding: analyze_results reports its leads as not laid out, and does not throw",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["core"]["functionalDescription"]["shape"] = "DR 5.3x3.6 + RIS 7.3x3.5x5.7";
    auto magnetic = autocomplete(fixture["magnetic"]);
    REQUIRE(magnetic.get_mutable_coil().is_real_winding_blocking_applied());
    REQUIRE(magnetic.get_mutable_coil().get_turns_description());
    // Laying the leads out is still refused...
    REQUIRE_THROWS_WITH(magnetic.get_mutable_coil().get_connection_layout(), ContainsSubstring("has no lead exit"));
    // ... and analyze_results says so in one finding instead of throwing it.
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    std::vector<ResultFinding> findings;
    REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(mas));
    size_t notLaidOut = 0;
    for (const auto& finding : findings) {
        CAPTURE(finding.message);
        CHECK(finding.code != ResultFindingCode::LEAD_EXIT_CLASHES);
        if (finding.code == ResultFindingCode::LEAD_ROUTE_NOT_MODELLED) {
            ++notLaidOut;
            CHECK_THAT(finding.message, ContainsSubstring("states no base plate (letter M)"));
            CHECK_THAT(finding.message, !ContainsSubstring("ABT"));
        }
    }
    CHECK(notLaidOut == 1);
    settings.reset();
}

// Alf, 2026-10-26: painting a catalogue drumRing (no base plate M) under real winding paints its core,
// bobbin and turns, and no terminal leads: the part states no lead exit, so its leads have no route.
TEST_CASE("A catalogue drumRing under real winding paints its core, bobbin and turns without leads, and does not throw",
          "[lead-exit][constructive-model][coil][painter][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["core"]["functionalDescription"]["shape"] = "DR 5.3x3.6 + RIS 7.3x3.5x5.7";
    auto magnetic = autocomplete(fixture["magnetic"]);
    REQUIRE(magnetic.get_mutable_coil().is_real_winding_blocking_applied());
    REQUIRE(magnetic.get_mutable_coil().get_turns_description());
    const size_t turns = magnetic.get_mutable_coil().get_turns_description()->size();
    auto outFile = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");
    std::filesystem::create_directories(outFile);
    outFile.append("Test_Lead_Exit_Catalogue_DrumRing_Painted.svg");
    std::filesystem::remove(outFile);
    Painter painter(outFile);
    REQUIRE_NOTHROW(painter.paint_magnetic(magnetic));
    painter.export_svg();
    REQUIRE(std::filesystem::exists(outFile));
    std::ifstream svgFile(outFile);
    const std::string svg((std::istreambuf_iterator<char>(svgFile)), std::istreambuf_iterator<char>());
    CHECK_THAT(svg, ContainsSubstring("ferrite"));
    CHECK_THAT(svg, ContainsSubstring("bobbin"));
    CHECK_THAT(svg, ContainsSubstring("copper"));
    for (size_t turn = 0; turn < turns; ++turn) {
        CAPTURE(turn);
        CHECK_THAT(svg, ContainsSubstring(magnetic.get_mutable_coil().get_turns_description().value()[turn].get_name()));
    }
    CHECK_THAT(svg, !ContainsSubstring("class=\"connection_terminal\""));
    settings.reset();
}

TEST_CASE("A closed core's lead exit refuses what it cannot build",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    SECTION("no pad data on the winding's connections") {
        auto fixture = load_fixture("lead_exit_drumring_7332.json");
        fixture["magnetic"]["coil"]["functionalDescription"][0].erase("connections");
        REQUIRE_THROWS_WITH(autocomplete(fixture["magnetic"]).get_mutable_coil().get_connection_layout(),
                            ContainsSubstring("padWidth and padDepth"));
    }
    SECTION("the lead frame is thicker than the molded body's base") {
        auto fixture = load_fixture("lead_exit_molded_4030.json");
        // Base (B - D)/2 = 0.1 mm, under the 0.127 mm sourced frame stock.
        fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"]["B"]["nominal"] = 0.0023;
        REQUIRE_THROWS_WITH(autocomplete(fixture["magnetic"]).get_mutable_coil().get_connection_layout(),
                            ContainsSubstring("thicker than the molded body's base"));
    }
    settings.reset();
}

namespace {

// The leads of a drum-ring part are drawn to their electrodes; returns how many terminal routes there are.
size_t check_leads_drawn_to_electrodes(const ConnectionLayout& layout) {
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        CHECK(route.exitRunEnd == ConnectionRoute::ExitRunEnd::ELECTRODE);
        CHECK(route.exitRun.size() >= 2);
        CHECK(route.exitRunArc);
    }
    return terminals;
}

size_t count_ring_bore_clashes(const ConnectionLayout& layout) {
    size_t pastBore = 0;
    for (const auto& clash : layout.leadExitClashes) {
        if (clash.kind == LeadExitClash::Kind::LANE_PAST_RING_BORE) {
            ++pastBore;
            CHECK_THAT(clash.description, ContainsSubstring("past the ring bore"));
        }
    }
    return pastBore;
}

}  // namespace

TEST_CASE("A drum-ring lead lane that passes the ring bore is drawn and reported",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    SECTION("the band between the winding and the ring bore is too narrow for the lane") {
        settings.set_coil_use_real_winding_geometry(false);
        auto fixture = load_fixture("lead_exit_drumring_7332.json");
        // 10 turns of 0.355 mm wire (0.384 mm pitch), ideal winding: two layers whose outer face is at
        // 2.017 mm. The lane sits at 2.592 mm (outer face + 0.384 mm ride-over + half its pitch) and
        // with its 0.192 mm bend reaches 2.784 mm, past a ring bore brought in to K/2 = 2.675 mm.
        auto& winding = fixture["magnetic"]["coil"]["functionalDescription"][0];
        winding["wire"] = "Round 0.355 - Grade 1";
        winding["numberTurns"] = 10;
        winding["numberParallels"] = 1;
        fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"]["K"]["nominal"] = 0.00535;
        auto magnetic = autocomplete(fixture["magnetic"]);
        ConnectionLayout layout;
        REQUIRE_NOTHROW(layout = magnetic.get_mutable_coil().get_connection_layout());
        CHECK(check_leads_drawn_to_electrodes(layout) == 2);
        CHECK(count_ring_bore_clashes(layout) == 2);
    }
    SECTION("a bifilar end's two lanes and bend do not fit inside the ring bore under real winding") {
        settings.set_coil_use_real_winding_geometry(true);
        auto fixture = load_fixture("lead_exit_drumring_7332.json");
        // 4 bifilar turns of 0.315 mm wire (0.341 mm pitch): real winding lays the wires in three
        // layers with an outer face at 2.274 mm. The outer lane of an end sits at 2.787 mm and with
        // its 0.328 mm bend reaches 3.115 mm, past the 2.85 mm ring bore: one clash per lead.
        auto& winding = fixture["magnetic"]["coil"]["functionalDescription"][0];
        winding["wire"] = "Round 0.315 - Grade 1";
        winding["numberTurns"] = 4;
        winding["numberParallels"] = 2;
        auto magnetic = autocomplete(fixture["magnetic"]);
        REQUIRE(magnetic.get_mutable_coil().is_real_winding_blocking_applied());
        ConnectionLayout layout;
        REQUIRE_NOTHROW(layout = magnetic.get_mutable_coil().get_connection_layout());
        CHECK(check_leads_drawn_to_electrodes(layout) == 4);
        CHECK(count_ring_bore_clashes(layout) == 4);
    }
    settings.reset();
}

TEST_CASE("A seated drumRing closes its bottom through the seat joint, not a radial clearance",
          "[lead-exit][constructive-model][core][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    OpenMagnetics::Core seated(fixture["magnetic"]["core"]);
    auto gapping = seated.get_functional_description().get_gapping();
    REQUIRE(gapping.size() == 2);
    // Top: the radial clearance (K - A)/2 across the top flange; bottom: the residual seat joint.
    const auto& top = gapping[0].get_coordinates().value()[1] > 0 ? gapping[0] : gapping[1];
    const auto& seat = gapping[0].get_coordinates().value()[1] > 0 ? gapping[1] : gapping[0];
    CHECK_THAT(top.get_length(), WithinAbs((dimension(fixture, "K") - dimension(fixture, "A")) / 2, 1e-12));
    CHECK_THAT(seat.get_length(), WithinAbs(Constants().residualGap, 1e-15));
    CHECK_THAT(seat.get_area().value(), WithinAbs(std::numbers::pi / 4 * (std::pow(dimension(fixture, "J"), 2) - std::pow(dimension(fixture, "K"), 2)), 1e-9));
    // The envelope is the square plate.
    CHECK_THAT(seat.get_coordinates().value()[1], WithinAbs(-dimension(fixture, "B") / 2 + dimension(fixture, "F"), 1e-12));
    auto processed = seated.get_processed_description().value();
    CHECK_THAT(processed.get_width(), WithinAbs(dimension(fixture, "M"), 1e-12));
    CHECK_THAT(processed.get_depth(), WithinAbs(dimension(fixture, "M"), 1e-12));
    // Magnetically: the seat joint (Constants::residualGap across the ring annulus) replaces the
    // bottom flange's radial clearance, so the seated part's inductance is the higher one. Like for
    // like: the same drum and ring letters in the catalogue construction (no plate; the ring spans
    // both flanges, L = B).
    auto seatedMagnetic = autocomplete(fixture["magnetic"]);
    auto catalogueFixture = fixture;
    auto& catalogueDimensions = catalogueFixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"];
    catalogueDimensions.erase("M");
    catalogueDimensions.erase("N");
    catalogueDimensions["L"]["nominal"] = dimension(fixture, "B");
    auto catalogueMagnetic = autocomplete(catalogueFixture["magnetic"]);
    const double seatedInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(seatedMagnetic)
                                        .get_magnetizing_inductance().get_nominal().value();
    const double catalogueInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(catalogueMagnetic)
                                           .get_magnetizing_inductance().get_nominal().value();
    CHECK(std::isfinite(seatedInductance));
    CHECK(seatedInductance > catalogueInductance);
    // A plate that does not carry its ring is refused.
    auto small = fixture["magnetic"]["core"];
    small["functionalDescription"]["shape"]["dimensions"]["M"]["nominal"] = 0.0070;
    REQUIRE_THROWS_WITH(OpenMagnetics::Core(small), ContainsSubstring("must carry the ring"));
    // Chamfers that leave no corner outside the ring are refused.
    auto chamfered = fixture["magnetic"]["core"];
    chamfered["functionalDescription"]["shape"]["dimensions"]["N"]["nominal"] = 0.0022;
    REQUIRE_THROWS_WITH(OpenMagnetics::Core(chamfered), ContainsSubstring("no corner area"));
}

// The lane round the winding clears MKF's OWN ride-over. On the seated drumRing 7332 fixture a return
// is laid under the outer layer (a Z dragback: 7 + 2 turns), so the outer face as wound stands one
// ride (one OD) further out than the layers' outer face; a consumer lifting the attached turn's end
// by that ride must still find the lead one coated diameter outside it. Before the lane took the
// ride, exitRun[0] sat ~3 um outside the lifted attach.
TEST_CASE("A closed core's lead lane clears the ride-over MKF's own returns impose on the winding's outer face",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::RING_NOTCH);
    const auto layout = coil.get_connection_layout();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    double windingOuterRadius = std::numeric_limits<double>::lowest();
    const auto layers = coil.get_layers_description().value();
    for (const auto& layer : layers) {
        windingOuterRadius = std::max(windingOuterRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    }
    const auto turns = coil.get_turns_description().value();
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        CAPTURE(entrance, route.side);
        const double ride = layout.ride_at(windingOuterRadius, route.side);
        CAPTURE(ride, windingOuterRadius, diameter);
        REQUIRE(ride > 0);   // anti-vacuous: this fixture's lanes DO have a ride to clear
        REQUIRE(route.exitRunRide);
        CHECK_THAT(route.exitRunRide.value(), WithinAbs(ride, 1e-12));
        REQUIRE(route.exitRunArc);
        CHECK(route.exitRunArc->radius >= windingOuterRadius + ride + diameter / 2 - 1e-9);
        const std::string turnName = entrance ? route.toTurn : route.fromTurn;
        std::optional<double> attachedRadius;
        for (const auto& turn : turns) {
            if (turn.get_name() == turnName) {
                attachedRadius = turn.get_coordinates()[0];
            }
        }
        REQUIRE(attachedRadius);
        REQUIRE(route.exitRun.size() >= 1);
        const double startRadius = std::hypot(route.exitRun[0][0], route.exitRun[0][2]);
        CAPTURE(attachedRadius.value(), startRadius);
        CHECK(startRadius - (attachedRadius.value() + ride) >= diameter - 1e-9);
    }
    CHECK(terminals == 2);
    settings.reset();
}

TEST_CASE("A molded body whose turns share one layer has no ride for its lead lane to clear",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_molded_4030.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::LEAD_FRAME);
    const auto layout = coil.get_connection_layout();
    // No return lays a lane (no Z dragback, no edge continuation), so nothing rides over anything.
    CHECK(layout.rideLevels.empty());
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        CHECK(route.kind != ConnectionKind::Z_DRAGBACK);
        CHECK(route.kind != ConnectionKind::EDGE_CONTINUATION);
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        REQUIRE(route.exitRunRide);
        CHECK(route.exitRunRide.value() == 0);
    }
    CHECK(terminals == 2);
    settings.reset();
}

// ABT #1424 / #1383: a bifilar Z winding's layer changes are DRAGBACKS -- each parallel finishes
// layer 0 at the top of the window and restarts layer 1 at the bottom -- and they are judged as
// what they are. The #1424 plane-link rule used to classify them by the slope between the two
// stations (dy <= dx + 2 OD K), which admitted this 1.3 mm return over a 0.27 mm radial step as a
// "link" and refused the two parallels 0.0577 mm apart in the connection plane, where nothing is
// drawn: a dragback steps out radially at the height its layer finished and runs back along the
// lane at the next layer's radius, and sibling dragbacks lie side by side along that lane.
//
// Everything checked here is re-derived from the stations, not read from the rule it guards:
//   - the return is a dragback because the conductor climbed through layer 0 and goes BACK the
//     other way to start layer 1 -- the definition, from the turn order alone;
//   - the straight station-to-station diagonals of the two parallels pass closer than one OD (the
//     number the slope-classified rule refused), while the drawn chains share the lane radius and
//     overlap along it, so the two need lanes of their own around the column;
//   - two lanes around a round column need 2 asin(OD / r) of azimuth at the innermost radius r
//     the chains occupy, and a revolution has 2 pi;
//   - MKF's ride model agrees they lie side by side: ONE ride level, one OD tall, at the lane.
TEST_CASE("A bifilar drum-ring's Z dragbacks are laid side by side along their lane, not judged as plane links",
          "[lead-exit][constructive-model][coil][abt1424-dragback][abt-1424]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.resolve_bobbin().get_processed_description().value().get_column_shape() == ColumnShape::ROUND);
    const auto turns = coil.get_turns_description().value();
    const double od = coil.get_wires()[0].get_maximum_outer_height();
    const double width = coil.get_wires()[0].get_maximum_outer_width();

    struct Return {
        std::string from, to;
        double xa, ya, xb, yb;
    };
    std::vector<Return> returns;
    for (int64_t parallel = 0; parallel < 2; ++parallel) {
        std::vector<const Turn*> sequence;
        for (const auto& turn : turns) {
            if (turn.get_parallel() == parallel) {
                sequence.push_back(&turn);
            }
        }
        for (size_t k = 1; k + 1 < sequence.size(); ++k) {
            const Turn& previous = *sequence[k - 1];
            const Turn& a = *sequence[k];
            const Turn& b = *sequence[k + 1];
            if (previous.get_layer() != a.get_layer() || a.get_layer() == b.get_layer()) {
                continue;
            }
            const double climb = a.get_coordinates()[1] - previous.get_coordinates()[1];
            const double back = b.get_coordinates()[1] - a.get_coordinates()[1];
            if (climb * back < 0) {
                returns.push_back({a.get_name(), b.get_name(), a.get_coordinates()[0], a.get_coordinates()[1],
                                   b.get_coordinates()[0], b.get_coordinates()[1]});
            }
        }
    }
    REQUIRE(returns.size() == 2);

    const auto layout = coil.get_connection_layout();
    for (const auto& r : returns) {
        size_t found = 0;
        for (const auto& route : layout.routes) {
            if (route.fromTurn == r.from && route.toTurn == r.to) {
                ++found;
                CHECK(route.kind == ConnectionKind::Z_DRAGBACK);
            }
        }
        CAPTURE(r.from, r.to);
        CHECK(found == 1);
    }

    const auto& p = returns[0];
    const auto& q = returns[1];
    const auto pointToSegment = [](double px, double py, double ax, double ay, double bx, double by) {
        const double vx = bx - ax, vy = by - ay;
        const double u = std::clamp(((px - ax) * vx + (py - ay) * vy) / (vx * vx + vy * vy), 0.0, 1.0);
        return std::hypot(px - (ax + u * vx), py - (ay + u * vy));
    };
    const double diagonals = std::min({pointToSegment(p.xa, p.ya, q.xa, q.ya, q.xb, q.yb),
                                       pointToSegment(p.xb, p.yb, q.xa, q.ya, q.xb, q.yb),
                                       pointToSegment(q.xa, q.ya, p.xa, p.ya, p.xb, p.yb),
                                       pointToSegment(q.xb, q.yb, p.xa, p.ya, p.xb, p.yb)});
    CAPTURE(diagonals, od);
    CHECK(diagonals < od);
    // The descents share the lane radius and overlap along it: one plane cannot hold both.
    CHECK_THAT(p.xb, WithinAbs(q.xb, 1e-9));
    CHECK(std::max(std::min(p.ya, p.yb), std::min(q.ya, q.yb)) < std::min(std::max(p.ya, p.yb), std::max(q.ya, q.yb)));
    const double innermost = std::min({p.xa, p.xb, q.xa, q.xb});
    CAPTURE(innermost);
    REQUIRE(od <= innermost);
    CHECK(2 * std::asin(od / innermost) <= 2 * std::numbers::pi);

    size_t levelsAtLane = 0;
    for (const auto& level : layout.rideLevels) {
        if (std::abs(level.radius - p.xb) <= 0.5 * width) {
            ++levelsAtLane;
            CHECK_THAT(level.height, WithinAbs(width, 1e-12));
        }
    }
    CHECK(levelsAtLane == 1);

    // ABT #2274 (Alf, 2026-11-05: "MKF owns the placement"): MKF lays the two parallels' returns side
    // by side itself -- parallel 0 on the connection plane, parallel 1 the one coated envelope off it
    // that two lines parallel to the axis need at the lane radius, asin(od / r) -- and the ride over
    // the level is the exact circle reach at that azimuth, not one OD.
    std::vector<double> azimuths;
    for (const auto& r : returns) {
        for (const auto& route : layout.routes) {
            if (route.fromTurn == r.from && route.toTurn == r.to) {
                REQUIRE(route.laneAzimuth.has_value());
                azimuths.push_back(route.laneAzimuth.value());
            }
        }
    }
    REQUIRE(azimuths.size() == 2);
    const double lane = p.xb;
    const double expectedAzimuth = std::asin(width / std::max(lane, width));
    CHECK_THAT(std::min(azimuths[0], azimuths[1]), WithinAbs(0.0, 1e-12));
    CHECK_THAT(std::max(azimuths[0], azimuths[1]), WithinAbs(expectedAzimuth, 1e-11));
    for (const auto& level : layout.rideLevels) {
        if (std::abs(level.radius - lane) <= 0.5 * width) {
            CHECK_THAT(level.offPlane, WithinAbs(expectedAzimuth, 1e-11));
            const double s = std::sin(level.offPlane);
            const double reach = level.radius * std::cos(level.offPlane) -
                                 std::sqrt((level.radius - width) * (level.radius - width) - level.radius * level.radius * s * s);
            CHECK_THAT(level.rise, WithinAbs(reach, 1e-15));
            CHECK(level.rise > width);
            CHECK_THAT(layout.ride_at(level.radius, level.side), WithinAbs(level.rise, 1e-15));
        }
    }
    settings.reset();
}

// The dragback rule refuses a lane that cannot hold its siblings. The same wound bifilar layout,
// moved radially inward (every station by the same amount, so the layout's kinds -- decided by
// the winding order and the layers -- are unchanged) until the source layer sits at 0.99 OD from
// the axis: two radial run-outs there cannot be separated by one diameter at any azimuth. (Not
// further in: the terminal lead leaving that layer must still ramp off its turn, which needs the
// turn's radius above the lead's planned bend radius, ABT #1336 -- and that refusal would come
// first.)
TEST_CASE("17_cllc: on a rectangular column MKF lays each conductor's return lane along the face, the primary group on -x and the other groups on +x (spec S0)",
          "[constructive-model][coil][real-winding][abt2274-lane-x]") {
    // ABT #2274 (Alf, 2026-11-06: "MKF owns each return lane's x along the face"): the x of a Z
    // dragback's lane on a flat face is MKF's (ConnectionRoute::laneOffset), not the 3D builder's. The
    // design: MAS example 17_cllc_xfmr_e5528_3c92a, E55/28 (rectangular column), Primary and Secondary
    // each wound in two layers with one Z dragback, no bobbin pins -- one face under S0.
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("..")
                           .append("MAS")
                           .append("examples")
                           .append("17_cllc_xfmr_e5528_3c92a.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coilIn);
    auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
    const auto bobbin = coil.resolve_bobbin();
    REQUIRE(bobbin.get_processed_description());
    REQUIRE(bobbin.get_processed_description()->get_column_shape() == MAS::ColumnShape::RECTANGULAR);
    const auto layout = coil.get_connection_layout();
    std::map<std::string, std::set<double>> lanesByWinding;
    double widest = 0;
    for (const auto& route : layout.routes) {
        INFO(label(route));
        CHECK(!route.laneAzimuth.has_value());   // a flat face has no azimuth
        if (route.kind != ConnectionKind::Z_DRAGBACK) {
            CHECK(!route.laneOffset.has_value());
            continue;
        }
        CHECK(route.side == 0);   // S0: one face
        REQUIRE(route.laneOffset.has_value());
        lanesByWinding[route.winding].insert(route.laneOffset.value());
        const auto wire = coil.resolve_wire(coil.get_winding_index_by_name(route.winding));
        widest = std::max(widest, wire.get_maximum_outer_width());
    }
    REQUIRE(lanesByWinding.size() == 2);
    REQUIRE(lanesByWinding.at("Primary").size() == 1);     // one lane per conductor
    REQUIRE(lanesByWinding.at("Secondary").size() == 1);
    const double primaryX = *lanesByWinding.at("Primary").begin();
    const double secondaryX = *lanesByWinding.at("Secondary").begin();
    CHECK(primaryX <= 0);      // the primary group on the -x half (lane 0 at the crossing)
    CHECK(secondaryX > 0);     // the other group on the +x half
    CHECK_THAT(secondaryX - primaryX, WithinAbs(widest, 1e-12));   // insulation touching
    settings.reset();
}

TEST_CASE("13_current_sense: the primary group's leads are laid first, so a Secondary listed before it keeps the face (spec S0, ABT #2274)",
          "[constructive-model][coil][real-winding][abt2274-s0-order]") {
    // ABT #2274 (Alf, 2026-11-06: "Primary first, move only misfits"). The design: MAS example
    // 13_current_sense_er95_n87, ER 9.5/5, no bobbin pins; the MAS lists the Secondary (80 turns,
    // 0.1 mm) BEFORE the one-turn Primary (0.5 mm). Laid in MAS order, the Secondary's entrance took
    // the crossing (x = 0), the Primary's one-turn legs, centred on the crossing, met it, and the
    // section-23 move put the Secondary -- whose leads fitted -- on the other face.
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("..")
                           .append("MAS")
                           .append("examples")
                           .append("13_current_sense_er95_n87.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
    REQUIRE(coilIn.get_functional_description().front().get_isolation_side() != MAS::IsolationSide::PRIMARY);   // the premise
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coilIn);
    auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto bobbin = coil.resolve_bobbin();
    REQUIRE(bobbin.get_processed_description());
    REQUIRE((!bobbin.get_processed_description()->get_pins() || bobbin.get_processed_description()->get_pins()->empty()));
    const auto layout = coil.get_connection_layout();
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        INFO(label(route));
        CHECK(route.side == 0);   // S0: ONE FACE, the Secondary included
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const auto& winding = coil.get_functional_description()[coil.get_winding_index_by_name(route.winding)];
        REQUIRE(route.laneHalf.has_value());
        CHECK(route.laneHalf.value() == (winding.get_isolation_side() == MAS::IsolationSide::PRIMARY ? -1 : 1));
        for (const auto& finding : route.findings) {
            CHECK(finding.kind != LeadFindingKind::LANES_DO_NOT_FIT);   // the face is not full
        }
    }
    CHECK(terminals == 4);
    settings.reset();
}

// ABT #2274 (Alf, 2026-11-07: "MKF checks height at the lane"). Under real winding a terminal lead on the
// FAR half of its turn's crossing -- an entrance upstream (x > 0), an exit past the crossing (x < 0) --
// puts its attach turn where the blocking never put it: by the height model every 3D builder draws
// (station - advance * along / L) the first turn stands below its station up to its lane, the last turn
// climbs above it to its lane. Measured here from the turns, the bobbin window and the routes alone: every
// such turn keeps its copper inside the winding window at its lane, and passes no other terminal lead's
// in-window run closer than the two coated sections -- or its lead carries LANES_DO_NOT_FIT. The designs
// are the UNMODIFIED MAS examples where the gaps were found (MVB++ drawing MKF's exitSlot, agent 9):
// 02_flyback_efd25_3c95 and 17_cllc_xfmr_e5528_3c92a (Primary exit past the crossing climbed 18 um /
// 46.6 um into the bobbin flange), 04_forward_xfmr_e3216_n87 (the 3.3V p0 entrance at +1.876 mm crossed
// the 12V entrance's run 18.6 um into it).
TEST_CASE("Real winding: a far-half terminal turn keeps its copper inside the window and off the leads it crosses, or its lead is reported (ABT #2274)",
          "[constructive-model][coil][real-winding][abt2274-lane-height]") {
    const std::string design = GENERATE(std::string("02_flyback_efd25_3c95"), std::string("17_cllc_xfmr_e5528_3c92a"),
                                        std::string("04_forward_xfmr_e3216_n87"));
    INFO(design);
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("..")
                           .append("MAS")
                           .append("examples")
                           .append(design + ".json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coilIn);
    auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto bobbin = coil.resolve_bobbin();
    REQUIRE(bobbin.get_processed_description());
    const auto window = bobbin.get_processed_description()->get_winding_windows().at(0);
    const double windowLow = window.get_coordinates().value().at(1) - window.get_height().value() / 2;
    const double windowHigh = window.get_coordinates().value().at(1) + window.get_height().value() / 2;
    const auto turns = coil.get_turns_description().value();
    const auto layout = coil.get_connection_layout();
    auto turn_index = [&](const std::string& name) {
        for (size_t t = 0; t < turns.size(); ++t) {
            if (turns[t].get_name() == name) {
                return t;
            }
        }
        FAIL("turn '" << name << "' not found");
        return size_t{0};
    };
    struct Lead {
        const OpenMagnetics::ConnectionRoute* route;
        OpenMagnetics::TerminalTurnPath path;
        double slot;
        double along;
        double station;
        double advance;
        double radius;
        bool reported;
    };
    std::vector<Lead> leads;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route) || !route.exitSlot || route.axialLead) {
            continue;
        }
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        const size_t t = turn_index(entrance ? route.toTurn : route.fromTurn);
        Lead lead{&route, coil.terminal_turn_path(route, turns[t], layout.routes), route.exitSlot.value(), 0, turns[t].get_coordinates().at(1),
                  0, (entrance ? route.waypoints.back() : route.waypoints.front()).at(0), false};
        // The advance in the wire's own direction (to the next turn of the conductor; from the previous one for an exit).
        auto same_conductor = [&](size_t n) {
            return turns[n].get_winding() == turns[t].get_winding() && turns[n].get_parallel() == turns[t].get_parallel();
        };
        if (entrance) {
            for (size_t n = t + 1; n < turns.size(); ++n) {
                if (same_conductor(n)) {
                    lead.advance = turns[n].get_coordinates().at(1) - lead.station;
                    break;
                }
            }
        }
        else {
            for (size_t n = t; n-- > 0;) {
                if (same_conductor(n)) {
                    lead.advance = lead.station - turns[n].get_coordinates().at(1);
                    break;
                }
            }
        }
        lead.along = OpenMagnetics::Coil::terminal_lead_along(route, lead.path, lead.slot);
        for (const auto& finding : route.findings) {
            lead.reported = lead.reported || finding.kind == LeadFindingKind::LANES_DO_NOT_FIT;
        }
        leads.push_back(lead);
    }
    REQUIRE(leads.size() >= 4);
    size_t farHalf = 0;
    for (const auto& lead : leads) {
        const bool entrance = lead.route->kind == ConnectionKind::TERMINAL_ENTRANCE;
        const bool far = entrance ? lead.along > 0 : lead.along < 0;
        farHalf += far ? 1 : 0;
        INFO(label(*lead.route) << " lane x = " << lead.slot * 1e3 << " mm, station " << lead.station * 1e3 << " mm");
        const double half = (lead.path.turnRectangular ? lead.path.turnAxialExtent : lead.path.turnWireDiameter) / 2;
        const double radialHalf = (lead.path.turnRectangular ? lead.path.turnRadialExtent : lead.path.turnWireDiameter) / 2;
        auto height = [&](double a) { return lead.station - lead.advance * a / lead.path.length; };
        const double atLane = height(lead.along);
        const bool inWindow = !far || atLane + half <= std::max(windowHigh, lead.station + half) + 1e-9 &&
                              atLane - half >= std::min(windowLow, lead.station - half) - 1e-9;
        bool offLeads = true;
        // The turn's copper on the face: an entrance's first turn from its lane on, an exit's last turn up to it.
        const double faceAlong = OpenMagnetics::Coil::terminal_lane_along(lead.path, lead.path.faceHalfStraight + lead.path.cornerRadius);
        const double lo = entrance ? -faceAlong : lead.along;
        const double hi = entrance ? lead.along : faceAlong;
        for (const auto& other : leads) {
            if (other.route == lead.route || other.route->side != lead.route->side) {
                continue;
            }
            const double ext = other.path.stubRectangular ? std::max(other.path.stubAxialExtent, other.path.stubRadialExtent) : other.path.stubDiameter;
            const double aFrom = std::max(OpenMagnetics::Coil::terminal_lane_along(lead.path, other.slot - ext / 2), lo);
            const double aTo = std::min(OpenMagnetics::Coil::terminal_lane_along(lead.path, other.slot + ext / 2), hi);
            if (aFrom > aTo || std::abs(other.slot) < ext / 2) {
                continue;   // outside the turn's copper, or met by the blocking at the crossing
            }
            const double rHalf = (other.path.stubRectangular ? other.path.stubRadialExtent : other.path.stubDiameter) / 2;
            const double yHalf = (other.path.stubRectangular ? other.path.stubAxialExtent : other.path.stubDiameter) / 2;
            std::optional<double> atStation;
            for (const double a : {0.0, aFrom, aTo}) {
                const double y = height(a);
                double deepest = 0;
                for (size_t w = 0; w + 1 < other.route->waypoints.size(); ++w) {
                    const auto& p = other.route->waypoints[w];
                    const auto& q = other.route->waypoints[w + 1];
                    const double radialOverlap = std::min(lead.radius + radialHalf, std::max(p.at(0), q.at(0)) + rHalf) -
                                                 std::max(lead.radius - radialHalf, std::min(p.at(0), q.at(0)) - rHalf);
                    const double axialOverlap = std::min(y + half, std::max(p.at(1), q.at(1)) + yHalf) -
                                                std::max(y - half, std::min(p.at(1), q.at(1)) - yHalf);
                    if (radialOverlap > 0 && axialOverlap > 0) {
                        deepest = std::max(deepest, std::min(radialOverlap, axialOverlap));
                    }
                }
                if (!atStation) {
                    atStation = deepest;   // what the blocking holds at the crossing
                    continue;
                }
                if (deepest - atStation.value() > 1e-9) {
                    UNSCOPED_INFO("crosses " << label(*other.route) << " at y = " << y * 1e3 << " mm by " << (deepest - atStation.value()) * 1e6 << " um");
                    offLeads = false;
                }
            }
        }
        CHECK(((inWindow && offLeads) || lead.reported));
    }
    WARN(design << ": " << leads.size() << " terminal leads judged, " << farHalf << " on the far half of their crossing");
    // Every lane stays on its face (no lane beyond the face's straight).
    for (const auto& lead : leads) {
        INFO(label(*lead.route));
        CHECK(std::abs(lead.slot) <= lead.path.faceHalfStraight + lead.path.cornerRadius + 1e-12);
    }
    settings.reset();
}

TEST_CASE("two_switch_forward: the foil Secondary's axial leads, judged at their real extent, do not fit the primary's face, so the Secondary moves (spec S0, ABT #2274 / #2291)",
          "[constructive-model][coil][real-winding][abt2274-foil-half]") {
    // ABT #2274 / #2291 (Alf, 2026-11-06). The design: MAS example two_switch_forward_transformer_complete,
    // E 42/21/15 (rectangular column), Primary 80 turns of litz, Secondary 2 turns x 8 parallels of Foil 0.2
    // whose ends leave AXIALLY on a stated Round 0.5 lead (ABT #1907). Each end's 8 leads stand side by side
    // from the crossing, one 0.534 mm lead each: the starts over x = -4.272 .. 0 mm, the finishes over
    // 0 .. +4.272 mm, both halves of a face whose straight is 5.975 mm. The axial leads were left out of the
    // S0 judgement: they got no half, no LANES_DO_NOT_FIT, and stayed on the Primary's face, their starts on
    // the primary group's half.
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("..")
                           .append("MAS")
                           .append("examples")
                           .append("complete")
                           .append("two_switch_forward_transformer_complete.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coilIn);
    auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto layout = coil.get_connection_layout();
    std::set<int> primaryFaces, secondaryFaces;
    size_t axialLeads = 0;
    for (const auto& route : layout.routes) {
        INFO(label(route));
        if (!is_terminal(route)) {
            continue;
        }
        const auto& winding = coil.get_functional_description()[coil.get_winding_index_by_name(route.winding)];
        const bool primary = winding.get_isolation_side() == MAS::IsolationSide::PRIMARY;
        (primary ? primaryFaces : secondaryFaces).insert(route.side);
        for (const auto& finding : route.findings) {
            CHECK(finding.kind != LeadFindingKind::LANES_DO_NOT_FIT);   // settled: every group fits its face
        }
        if (primary) {
            REQUIRE_FALSE(route.axialLead);
            continue;
        }
        REQUIRE(route.axialLead);   // the premise: the foil's ends leave axially (ABT #1907)
        ++axialLeads;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        // Alone on its face (spec section 24): starts on x <= 0, finishes on x >= 0, where the sheets put them.
        REQUIRE(route.laneHalf.has_value());
        CHECK(route.laneHalf.value() == (entrance ? -1 : 1));
        CHECK(route.laneOverflow == std::optional<bool>(false));
        // ABT #2291: on the window edge's row the lead turns 90 degrees out of the face, at the route's planned
        // bend; its copper stays inside the window -- it never crosses the flange (Alf, 2026-11-06).
        const auto& axial = route.axialLead.value();
        const double outwardEdge = axial.windowEdge > axial.attachTo ? 1.0 : -1.0;
        CHECK(outwardEdge * (axial.windowEdge - axial.end) >= axial.outerDiameter / 2 - 1e-12);
        REQUIRE(route.waypoints.size() == 3);
        const auto& corner = route.waypoints[1];
        const auto& out = entrance ? route.waypoints.front() : route.waypoints.back();
        const auto& attach = entrance ? route.waypoints.back() : route.waypoints.front();
        CHECK_THAT(attach[0], Catch::Matchers::WithinAbs(axial.layerAxis, 1e-12));
        CHECK_THAT(attach[1], Catch::Matchers::WithinAbs(axial.attachFrom, 1e-12));
        CHECK_THAT(corner[0], Catch::Matchers::WithinAbs(axial.layerAxis, 1e-12));
        CHECK_THAT(corner[1], Catch::Matchers::WithinAbs(axial.end, 1e-12));
        CHECK_THAT(out[0], Catch::Matchers::WithinAbs(axial.radialEnd, 1e-12));
        CHECK_THAT(out[1], Catch::Matchers::WithinAbs(axial.end, 1e-12));
        CHECK(axial.radialEnd > axial.layerAxis + route.plannedBendRadius);
        CHECK(route.plannedBendRadius > 0);
        CHECK_THAT(route.routedLength, Catch::Matchers::WithinAbs(axial.length + axial.radialLength, 1e-9));
    }
    CHECK(axialLeads == 16);
    REQUIRE(primaryFaces.size() == 1);
    REQUIRE(secondaryFaces.size() == 1);
    CHECK(*primaryFaces.begin() == 0);   // the face's first group keeps it
    CHECK(*secondaryFaces.begin() == 1);  // the Secondary, whose leads do not fit its half there, moves
    settings.reset();
}

TEST_CASE("Sibling Z dragbacks whose lane cannot hold them side by side are refused",
          "[constructive-model][coil][abt1424-dragback][abt-1424]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    const double od = coil.get_wires()[0].get_maximum_outer_height();
    // As wound: accepted.
    REQUIRE_NOTHROW(coil.refuse_sibling_layer_links_closer_than_the_wire());
    auto turns = coil.get_turns_description().value();
    double innermost = std::numeric_limits<double>::max();
    for (const auto& turn : turns) {
        innermost = std::min(innermost, turn.get_coordinates()[0]);
    }
    const double shift = innermost - 0.99 * od;
    REQUIRE(shift > 0);
    for (auto& turn : turns) {
        auto coordinates = turn.get_coordinates();
        coordinates[0] -= shift;
        turn.set_coordinates(coordinates);
    }
    coil.set_turns_description(turns);
    REQUIRE_THROWS_WITH(coil.refuse_sibling_layer_links_closer_than_the_wire(),
                        ContainsSubstring("Z dragbacks") && ContainsSubstring("side by side along the lane") &&
                            ContainsSubstring("no azimuth separates two radial run-outs by one diameter"));
    settings.reset();
}

// ABT #1424 follow-up: the LINK-PITCH pass asks MKF's connection layout what each layer change is
// while the coil is still being wound. That question must not plan the terminal leads: on this EFD
// 20/10/7 (IRREGULAR column) the half-wound layout has a Secondary terminal stub too short for its
// bends, and the ramp it would need has no surface curvature to be planned on -- an error that
// belongs to the finished coil, whose own leads lay out fine here.
TEST_CASE("Asking what a layer change is while winding plans no terminal lead",
          "[constructive-model][coil][abt-1424]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("test_circuitsimulatorexporter_simba_json_87.json");
    std::optional<OpenMagnetics::Magnetic> magnetic;
    // The fixture IS a magnetic (core + coil), not a MAS document.
    REQUIRE(fixture.contains("coil"));
    REQUIRE_NOTHROW(magnetic.emplace(autocomplete(fixture)));
    auto coil = magnetic->get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    std::vector<ConnectionRoute> routes;
    REQUIRE_NOTHROW(coil.get_connection_reserved_spaces(&routes));
    size_t terminalRoutes = 0;
    for (const auto& route : routes) {
        if (route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT) {
            ++terminalRoutes;
        }
    }
    // Every conductor's two ends are laid out on the finished coil: Primary 2 parallels, Secondary 1.
    CHECK(terminalRoutes == 6);
    settings.reset();
}


// Alf, 2026-10-31: a closed core's lead whose lane lies beyond its turn's half extent leaves the turn
// tangentially (ABT #1705) and keeps its radial run from the corner's end, so its exit run starts from
// that run at its slot, as any other lead's does.
TEST_CASE("A closed core's tangential terminal lead starts its exit run from its radial run",
          "[lead-exit][constructive-model][coil][abt-1336][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    auto magnetic = autocomplete(fixture["magnetic"]);
    const auto layout = magnetic.get_mutable_coil().get_connection_layout();
    size_t tangential = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route) || !route.tangentDeparture) {
            continue;
        }
        ++tangential;
        CHECK(route.tangentDeparture->kind != TangentDeparture::Kind::LATERAL_FACE);
        REQUIRE(route.exitSlot);
        REQUIRE_FALSE(route.exitRun.empty());
        // The exit run starts on the radial run, at the slot, beyond the corner's end.
        CHECK_THAT(route.exitRun.front()[0], WithinAbs(route.exitSlot.value(), 1e-12));
        CHECK(-route.exitRun.front()[2] > -route.tangentDeparture->cornerEnd[2]);
    }
    CHECK(tangential > 0);
    settings.reset();
}

// Alf, 2026-10-31: an IRREGULAR column states no surface curvature, but MKF winds it as a rounded
// rectangle (Coil::terminal_turn_path), so a terminal ramp off it is planned on that turn's own curvature
// where the lead leaves it, and the lanes of a coil wound on one are laid out.
TEST_CASE("A terminal ramp off an IRREGULAR column is planned on the turn MKF wound there",
          "[constructive-model][coil][abt-1336]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto magnetic = autocomplete(load_fixture("test_circuitsimulatorexporter_simba_json_87.json"));
    auto coil = magnetic.get_coil();
    REQUIRE(coil.resolve_bobbin().get_processed_description()->get_column_shape() == ColumnShape::IRREGULAR);
    std::vector<ConnectionRoute> routes;
    REQUIRE_NOTHROW(coil.get_connection_reserved_spaces(&routes));
    // The lane search plans ramps at its candidate slots (none is left at the slots it takes); before
    // this rule every such plan threw for the column's shape.
    size_t terminals = 0;
    for (const auto& route : routes) {
        terminals += is_terminal(route);
    }
    CHECK(terminals == 6);
    settings.reset();
}

}  // namespace lead_exit_suite

// ABT #1397: the seated drum-ring's coil ends are joined to their terminal legs by a BOND solid (the
// wire's own conductor, thermocompression-bonded): from the leg's top face (a shared face) up to the
// wire's centreline, so the bare conductor -- one coating thickness above the leg -- reaches the metal.
TEST_CASE("A seated drum-ring's coil ends are bonded to their terminal legs through the wire's coating gap",
          "[lead-exit][constructive-model][coil][abt-1397]") {
    using namespace lead_exit_suite;
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    auto fixture = lead_exit_suite::load_fixture("lead_exit_drumring_7332.json");
    auto magnetic = lead_exit_suite::autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    const auto layout = coil.get_connection_layout();
    const double coated = coil.get_wires()[0].get_maximum_outer_width();
    const double bare = coil.get_wires()[0].get_maximum_conducting_width();
    REQUIRE(bare < coated);
    size_t bonds = 0;
    for (const auto& route : layout.routes) {
        if (!lead_exit_suite::is_terminal(route)) {
            continue;
        }
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        const LeadExitSolid* bond = nullptr;
        const LeadExitSolid* leg = nullptr;
        for (const auto& solid : layout.leadExitSolids) {
            if (solid.entrance != entrance) continue;
            if (solid.kind == LeadExitSolid::Kind::BOND && solid.parallel == route.parallel) bond = &solid;
            if (solid.kind == LeadExitSolid::Kind::ELECTRODE && solid.name.find("terminal leg") != std::string::npos) leg = &solid;
        }
        REQUIRE(bond != nullptr);
        REQUIRE(leg != nullptr);
        ++bonds;
        // Shared face: the bond stands on the leg's top face.
        CHECK(std::abs(bond->body.y0 - leg->body.y1) < 1e-12);
        // It reaches the wire's centreline (the run's last point) and so overlaps the bare conductor,
        // whose bottom is one coating thickness above the leg.
        const auto& end = route.exitRun.back();
        CHECK(std::abs(bond->body.y1 - end[1]) < 1e-12);
        CHECK(bond->body.y1 > leg->body.y1 + (coated - bare) / 2);
        // Volume > 0 (shoelace area x height), and the wire's own conductor.
        double area = 0;
        const auto& xz = bond->body.xz;
        for (size_t k = 0; k < xz.size(); ++k) {
            const auto& p = xz[k];
            const auto& q = xz[(k + 1) % xz.size()];
            area += p[0] * q[1] - q[0] * p[1];
        }
        CHECK(std::abs(area) / 2 * (bond->body.y1 - bond->body.y0) > 0);
        CHECK(bond->material == coil.resolve_wire(0).resolve_material().get_name());
    }
    CHECK(bonds == 2);
    settings.reset();
}

TEST_CASE("terminal_exit_slots: tied members take the order whose turns clear each other's leads at the lanes they pass (ABT #1336)",
          "[constructive-model][coil][abt-1336]") {
    // 14_dab's Secondary exits (mm -> m): four parallels leave from their own rows (no stubs) on the outer
    // layer, r = 21.4375 mm, whose one turn per parallel is a FINAL_LANDING descending 37.557 mm per
    // revolution. Each last turn arrives at its lane from +x, passing the lanes beyond its own, and at a
    // lane x it stands 37.557 asin(x / r) / 2 pi above its row; the drawn lead leaves its turn there and
    // drops along the turn axis to its row. With the lowest row first, parallel 0's landing passes
    // parallel 1's lane 0.658 mm under that drop's foot, against the 0.855 mm coated envelope.
    const double od = 0.855e-3;
    const double bend = 0.8275e-3;
    const double radius = 21.4375e-3;
    const std::vector<double> rows = {-18.8410e-3, -17.9443e-3, -17.0477e-3, -16.1510e-3};
    std::vector<ConnectionRoute> routes;
    for (int64_t parallel = 0; parallel < 4; ++parallel) {
        ConnectionRoute r;
        r.winding = "Secondary";
        r.parallel = parallel;
        r.kind = ConnectionKind::TERMINAL_EXIT;
        r.side = 1;
        r.waypoints = {{radius, rows[parallel]}, {34.5025e-3, rows[parallel]}};
        r.plannedBendRadius = bend;
        r.fromTurn = "t" + std::to_string(parallel);
        // Spec section 21: each last turn arrives at its lane from +x, i.e. travels towards -x: the
        // clockwise turn of a window right of the column (Coil::terminal_wrap_direction).
        r.wrapDirection = -1;
        routes.push_back(r);
    }
    const std::vector<double> diameters(routes.size(), od);
    const WindingSides secondaryOnly = {{"Secondary", IsolationSide::PRIMARY}};
    auto slots_for = [&](const std::vector<double>& advance, bool realWinding) {
        return OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, rows, realWinding, advance,
                                                        realWinding ? round_paths(routes, diameters) : std::vector<std::optional<TerminalTurnPath>>{},
                                                        secondaryOnly);
    };
    const std::vector<double> descending(routes.size(), -37.557e-3);
    const std::vector<double> ascending(routes.size(), 3.6618e-3);
    SECTION("descending landings: the highest row stands first -- 14_dab's p3 / p2 / p1 / p0") {
        const auto slots = slots_for(descending, true);
        for (int64_t parallel = 0; parallel < 4; ++parallel) {
            CHECK(std::abs(slots[parallel].value() - (3 - parallel) * od) < 1e-12);
        }
    }
    SECTION("the clearance is measured, not signed: the ascending order fails by the landing's own height") {
        // What the order guards, computed here from the plan: parallel 0 at lane 0 passes parallel 1's lane od.
        const double dip = 37.557e-3 * std::asin(od / radius) / (2 * std::numbers::pi);
        const double found = (rows[1] - rows[0]) - dip;
        CHECK(found < od);   // 0.658 mm < 0.855 mm: the ascending order collides
        CHECK(OpenMagnetics::Coil::terminal_turn_height_at(routes[0], -37.557e-3, round_paths(routes, std::vector<double>(routes.size(), od))[0].value(), od) ==
              Catch::Approx(rows[0] + dip).margin(1e-15));
    }
    SECTION("ascending turns: the lowest row stands first, as under ideal winding") {
        const auto slots = slots_for(ascending, true);
        for (int64_t parallel = 0; parallel < 4; ++parallel) {
            CHECK(std::abs(slots[parallel].value() - parallel * od) < 1e-12);
        }
    }
    SECTION("ideal winding draws no turn past the crossing: ascending, whatever the advance") {
        const auto slots = slots_for(descending, false);
        for (int64_t parallel = 0; parallel < 4; ++parallel) {
            CHECK(std::abs(slots[parallel].value() - parallel * od) < 1e-12);
        }
    }
    SECTION("on a rectangular column the turn climbs over its whole path, so the highest row stands first too") {
        // A square column: 15 mm face straights either side of the crossing, 2 mm corners. MVB++ spreads
        // a rect turn's advance uniformly over its whole path from the crossing (appendRectWrap), so at
        // lane x an exit's last turn stands advance * x / length off its station -- not at the station.
        TerminalTurnPath rect{ColumnShape::RECTANGULAR, 8 * 15e-3 + 2 * std::numbers::pi * 2e-3, 15e-3, 2e-3, 17e-3, 13e-3};
        rect.turnWireDiameter = od;   // hand-built routes are unsleeved: turn and stub are the route's own wire
        rect.stubDiameter = od;
        const std::vector<std::optional<TerminalTurnPath>> paths(routes.size(), rect);
        CHECK(OpenMagnetics::Coil::terminal_turn_height_at(routes[0], -37.557e-3, rect, od) ==
              Catch::Approx(rows[0] + 37.557e-3 * od / rect.length).margin(1e-15));
        // Round the corner: the straight, then the arc's own length.
        const double cornerX = 15e-3 + 1e-3;
        CHECK(OpenMagnetics::Coil::terminal_turn_height_at(routes[0], -37.557e-3, rect, cornerX) ==
              Catch::Approx(rows[0] + 37.557e-3 * (15e-3 + 2e-3 * std::asin(0.5)) / rect.length).margin(1e-15));
        // An entrance wrap BEGINS at its slot past the crossing, on the -x half (ABT #1705): it stands
        // above its station there; upstream of the crossing it would stand below it, as an exit does.
        ConnectionRoute entrance = routes[0];
        entrance.kind = ConnectionKind::TERMINAL_ENTRANCE;
        entrance.toTurn = "t0";
        entrance.waypoints = {{34.5025e-3, rows[0]}, {radius, rows[0]}};
        CHECK(OpenMagnetics::Coil::terminal_turn_height_at(entrance, 3.6618e-3, rect, -od) ==
              Catch::Approx(rows[0] + 3.6618e-3 * od / rect.length).margin(1e-15));
        CHECK(OpenMagnetics::Coil::terminal_turn_height_at(entrance, 3.6618e-3, rect, od) ==
              Catch::Approx(rows[0] - 3.6618e-3 * od / rect.length).margin(1e-15));
        CHECK_THROWS_WITH(OpenMagnetics::Coil::terminal_turn_height_at(routes[0], -37.557e-3, rect, 17.5e-3),
                          Catch::Matchers::ContainsSubstring("beyond its turn's face"));
        const auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, rows, true, descending, paths, secondaryOnly);
        for (int64_t parallel = 0; parallel < 4; ++parallel) {
            CHECK(std::abs(slots[parallel].value() - (3 - parallel) * od) < 1e-12);
        }
    }
    SECTION("a lane beyond its turn leaves tangentially off the back of the turn, one corner onto its run (ABT #1705)") {
        // Round turn of radius r: the corner's centre lies on the mid plane one bend radius inside the
        // lane, so the straight is the direct common tangent of the turn and that corner circle.
        const auto path = round_paths(routes, std::vector<double>(routes.size(), od))[0].value();
        const double lane = radius + 2e-3;
        const auto plan = OpenMagnetics::Coil::plan_terminal_tangential_lead(routes[0], path, lane);
        const double centre = lane - bend;
        const double theta = std::acos((radius - bend) / centre);
        CHECK(plan.cornerAngle == Catch::Approx(theta).margin(1e-12));
        CHECK(plan.cornerAngle < std::numbers::pi / 2);
        CHECK(plan.straightLength == Catch::Approx(std::sqrt(centre * centre - (radius - bend) * (radius - bend))).margin(1e-12));
        CHECK(plan.along == Catch::Approx(radius * (std::numbers::pi / 2 + theta)).margin(1e-12));
        CHECK(plan.tangentPoint[1] > 0);   // off the BACK of the turn
        CHECK(std::hypot(plan.tangentPoint[0], plan.tangentPoint[1]) == Catch::Approx(radius).margin(1e-12));
        CHECK(plan.cornerCentre[0] == Catch::Approx(centre).margin(1e-15));
        // At exactly the turn's radius it leaves straight down the side: no straight, no corner, and
        // the same place on the turn as a lane on the face there.
        const auto side = OpenMagnetics::Coil::plan_terminal_tangential_lead(routes[0], path, radius);
        CHECK(side.cornerAngle == Catch::Approx(0.0).margin(1e-7));
        CHECK(side.straightLength == Catch::Approx(0.0).margin(1e-6));
        CHECK(side.along == Catch::Approx(OpenMagnetics::Coil::terminal_lane_along(path, radius)).margin(1e-9));
        // An entrance mirrors onto the -x half.
        ConnectionRoute entrance = routes[0];
        entrance.kind = ConnectionKind::TERMINAL_ENTRANCE;
        entrance.waypoints = {{34.5025e-3, rows[0]}, {radius, rows[0]}};
        const auto mirrored = OpenMagnetics::Coil::plan_terminal_tangential_lead(entrance, path, -lane);
        CHECK(mirrored.along == Catch::Approx(-plan.along).margin(1e-12));
        CHECK(mirrored.tangentPoint[0] == Catch::Approx(-plan.tangentPoint[0]).margin(1e-15));
        // Spec S0 (Alf, 2026-10-25, one face): an exit lane on the -x half is the FAR half for an exit,
        // and MKF plans its tangential lead too (it was refused as the "wrong half" before S0). On a
        // round turn, whose face has no straight, the turn runs on to the crossing and the lead leaves
        // it there, on the turn's centreline, level along the face line; one 90 degree corner of its
        // planned bend radius puts it on its radial run at the lane.
        OpenMagnetics::TerminalTangentialLead far;
        REQUIRE_NOTHROW(far = OpenMagnetics::Coil::plan_terminal_tangential_lead(routes[0], path, -lane));
        CHECK(far.farHalf);
        CHECK_FALSE(plan.farHalf);
        CHECK(far.lane == -lane);
        CHECK(far.tangentPoint[0] == Catch::Approx(0.0).margin(1e-15));
        CHECK(std::hypot(far.tangentPoint[0], far.tangentPoint[1]) == Catch::Approx(radius).margin(1e-12));
        CHECK(far.cornerAngle == Catch::Approx(std::numbers::pi / 2).margin(1e-15));
        CHECK(far.cornerEnd[0] == Catch::Approx(-lane).margin(1e-15));
        CHECK(std::hypot(far.cornerStart[0] - far.cornerCentre[0], far.cornerStart[1] - far.cornerCentre[1]) == Catch::Approx(bend).margin(1e-15));
        CHECK(std::hypot(far.cornerEnd[0] - far.cornerCentre[0], far.cornerEnd[1] - far.cornerCentre[1]) == Catch::Approx(bend).margin(1e-15));
        CHECK_THROWS_WITH(OpenMagnetics::Coil::plan_terminal_tangential_lead(routes[0], path, radius / 2),
                          Catch::Matchers::ContainsSubstring("on its turn's face"));
    }
    SECTION("tied members whose turns advance towards each other but stand apart cannot meet, so both take the crossing lane") {
        // Stacked sections: the rows stand 5 mm apart, and the lower turn rises while the upper falls. A
        // direction rule refuses this; the distance at the lanes passed is far beyond the envelope.
        // ABT #1705 (Alf, 2026-10-24, "lanes only where leads can meet", spec P5): leads that cannot meet
        // in 3D need no disjoint lanes, so both share the crossing lane. (Written for ascending lanes
        // 0 / od before that rule; the tie ORDER is guarded by the sections above, whose rows touch.)
        std::vector<double> farRows = {-20e-3, -15e-3};
        std::vector<ConnectionRoute> two = {routes[0], routes[1]};
        two[0].waypoints = {{radius, farRows[0]}, {34.5025e-3, farRows[0]}};
        two[1].waypoints = {{radius, farRows[1]}, {34.5025e-3, farRows[1]}};
        const auto slots = OpenMagnetics::Coil::terminal_exit_slots(two, {od, od}, farRows, true, {3.6618e-3, -3.6618e-3}, round_paths(two, {od, od}), secondaryOnly);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == 0.0);
    }
    SECTION("members whose stubs do not tie are ordered by their stubs; only the tie is permuted") {
        // Parallel 3 climbs 1 mm before its run, so it ties with nobody and keeps its place by its span.
        routes[3].waypoints = {{radius, rows[3]}, {radius, rows[3] + 1e-3}, {34.5025e-3, rows[3] + 1e-3}};
        std::vector<std::optional<double>> slots;
        REQUIRE_NOTHROW(slots = slots_for(descending, true));
        CHECK(slots[3].value() == 0.0);   // exits by span descending: the climbing one first
        CHECK(slots[2].value() < slots[1].value());
        CHECK(slots[1].value() < slots[0].value());
    }
    SECTION("when no order clears, the closest order is drawn and reported naming it, the pair and the shortfall") {
        // A level ring (advance 0) between descending landings: some landing always crosses its row.
        // Spec section 17 (Alf 2026-10-26): reported on every tied member (was CHECK_THROWS_WITH the same message).
        std::vector<double> advance = descending;
        advance[1] = 0.0;
        TerminalLanePlan plan;
        REQUIRE_NOTHROW(plan = OpenMagnetics::Coil::terminal_lane_plan(routes, diameters, rows, true, advance,
                                                                       round_paths(routes, diameters), secondaryOnly));
        REQUIRE(plan.findings.size() == routes.size());
        for (size_t index = 0; index < routes.size(); ++index) {
            CHECK(plan.slots[index].has_value());
            REQUIRE_FALSE(plan.findings[index].empty());
            CHECK(plan.findings[index].front().kind == LeadFindingKind::TIED_LEADS_COLLIDE);
            CHECK_THAT(plan.findings[index].front().witness,
                       Catch::Matchers::ContainsSubstring("No order of the tied exit leads of winding 'Secondary'") &&
                           Catch::Matchers::ContainsSubstring("mm short"));
        }
    }
}

TEST_CASE("a tie no order clears is drawn and reported and the plan returned (spec section 17, ABT #1336)",
          "[constructive-model][coil][abt-1705]") {
    // The "when no order clears" fixture above: 14_dab's Secondary exits, descending landings with a
    // level ring (advance 0) on parallel 1, which no order clears.
    const double od = 0.855e-3;
    const double bend = 0.8275e-3;
    const double radius = 21.4375e-3;
    const std::vector<double> rows = {-18.8410e-3, -17.9443e-3, -17.0477e-3, -16.1510e-3};
    std::vector<ConnectionRoute> routes;
    for (int64_t parallel = 0; parallel < 4; ++parallel) {
        ConnectionRoute r;
        r.winding = "Secondary";
        r.parallel = parallel;
        r.kind = ConnectionKind::TERMINAL_EXIT;
        r.side = 1;
        r.waypoints = {{radius, rows[parallel]}, {34.5025e-3, rows[parallel]}};
        r.plannedBendRadius = bend;
        r.fromTurn = "t" + std::to_string(parallel);
        routes.push_back(r);
    }
    const std::vector<double> diameters(routes.size(), od);
    const WindingSides secondaryOnly = {{"Secondary", IsolationSide::PRIMARY}};
    std::vector<double> advance(routes.size(), -37.557e-3);
    advance[1] = 0.0;
    // Spec section 17 (Alf 2026-10-26; was: setting off throws, setting on returns the slots and warns on
    // stderr): the plan is returned and every tied member carries the tie.
    TerminalLanePlan plan;
    REQUIRE_NOTHROW(plan = OpenMagnetics::Coil::terminal_lane_plan(routes, diameters, rows, true, advance,
                                                                   round_paths(routes, diameters), secondaryOnly));
    REQUIRE(plan.slots.size() == 4);
    for (int64_t parallel = 0; parallel < 4; ++parallel) {
        CHECK(plan.slots[parallel].has_value());
        REQUIRE_FALSE(plan.findings[parallel].empty());
        CHECK(plan.findings[parallel].front().kind == LeadFindingKind::TIED_LEADS_COLLIDE);
        CHECK_THAT(plan.findings[parallel].front().witness, Catch::Matchers::ContainsSubstring("No order of the tied exit leads"));
    }

    // Spec section 17 covers EVERY collision refusal (Alf 2026-10-26): each is drawn and reported on its lead.
    SECTION("plan_terminal_tangential_lead: a lane whose tangent meets the turn past its far side is planned anyway") {
        // A round turn of radius 2 mm and a lead bent at 3 mm: from lane x = 5 mm the straight tangent to the
        // back of the turn meets it at 120 deg past the side, 7.33 mm from the crossing, more than half the
        // turn's 12.57 mm.
        const double turnRadius = 2e-3;
        ConnectionRoute exit;
        exit.winding = "Secondary";
        exit.parallel = 0;
        exit.kind = ConnectionKind::TERMINAL_EXIT;
        exit.fromTurn = "t0";
        exit.waypoints = {{turnRadius, 0.0}, {turnRadius + 5e-3, 0.0}};
        exit.plannedBendRadius = 3e-3;
        TerminalTurnPath path{ColumnShape::ROUND, 2 * std::numbers::pi * turnRadius, 0.0, turnRadius, turnRadius};
        path.turnWireDiameter = 0.5e-3;
        path.stubDiameter = 0.5e-3;
        const double lane = 5e-3;
        // Was: setting off throws (and terminal_turn_height_along throws past half a turn), setting on plans
        // it and warns on stderr. Now planned, and the plan carries the witness the route reports.
        OpenMagnetics::TerminalTangentialLead plan;
        REQUIRE_NOTHROW(plan = OpenMagnetics::Coil::plan_terminal_tangential_lead(exit, path, lane));
        double height = 0;
        REQUIRE_NOTHROW(height = OpenMagnetics::Coil::terminal_turn_height_along(exit, 1e-3, path, plan.along));
        CHECK_THAT(plan.farSideWitness, Catch::Matchers::ContainsSubstring("past the far side of the turn"));
        CHECK(plan.crossesFarSide);
        CHECK(plan.along > path.length / 2);
        CHECK_THAT(plan.along, Catch::Matchers::WithinAbs(turnRadius * (std::numbers::pi / 2 + std::acos(-0.5)), 1e-12));
        CHECK_THAT(plan.cornerEnd[0], Catch::Matchers::WithinAbs(lane, 1e-15));
        CHECK(std::isfinite(height));
    }
}

TEST_CASE("14_dab: the Secondary's descending landings take their exit lanes highest row first (ABT #1336)",
          "[constructive-model][coil][real-winding][abt-1336]") {
    // Each Secondary parallel ends on a FINAL_LANDING that descends a whole window per revolution, so at
    // the lanes beyond its own a last turn stands well above its row. Ascending order put parallel 0's
    // landing 28 um inside parallel 1's coated envelope under MVB++'s certified gate; the order that
    // clears stands parallel 3 first. The order is the one MVB++ built and certified the part with;
    // the lanes stand 3.225 mm further out than that build's (6.977 mm for parallel 3): planned
    // pitch-true, the Secondary's parallel 2 ENTRANCE stands 1.624 mm off its row at that build's slot
    // (4.439 mm), short of the 1.655 mm its two bends need, so it ramps, and its ramp's span moves the
    // entrance bundle and the exit bundle after it out along the face.
    auto coil = wind_dab_fixture();
    const auto layout = coil.get_connection_layout();
    std::map<int64_t, double> exitSlot;
    for (const auto& route : layout.routes) {
        if (route.winding == "Secondary" && route.kind == ConnectionKind::TERMINAL_EXIT) {
            REQUIRE(route.exitSlot.has_value());
            exitSlot[route.parallel] = route.exitSlot.value();
        }
    }
    REQUIRE(exitSlot.size() == 4);
    // Spec S0 (Alf, 2026-10-25, ONE FACE; this bobbin has no pins, so section 12 does not apply): the
    // Secondary no longer has a face of its own. It shares face 0 with the Primary on the x >= 0 half,
    // so its entrance bundle starts one od (0.855 mm, C3: the coated radii) past the Primary's p3 exit
    // at the crossing instead of at 0.
    // Spec P5 (lanes only where leads can meet; Alf, 2026-10-28, "they look fine"): the planner's 3D envelopes
    // judge Secondary exit p3 unable to meet entrance p0 (its turns stand 2.3 mm higher on the face), so it
    // shares p0's lane at 0.855 mm and p2 stands one od beyond it; p1 and p0 stay past the entrance bundle
    // (entrances at 0.855 / 4.30861 / 8.47408 / 9.32908 mm). The order (p3 first) is unchanged.
    // ABT #2117 (32e9130b) anchors the helical end bundles at their outer turn's edge, a few um off
    // the centre they were laid at, so p1 and p0 stand 5.53 um further out than the 10.18408 /
    // 11.03908 mm planned before it (Alf, 2026-10-31).
    CHECK_THAT(exitSlot[3], Catch::Matchers::WithinAbs(0.855e-3, 0.5e-6));
    CHECK_THAT(exitSlot[2], Catch::Matchers::WithinAbs(1.71e-3, 0.5e-6));
    CHECK_THAT(exitSlot[1], Catch::Matchers::WithinAbs(10.18961e-3, 0.5e-6));
    CHECK_THAT(exitSlot[0], Catch::Matchers::WithinAbs(11.04461e-3, 0.5e-6));
    Settings::GetInstance().reset();
}

// ABT #1705 (Alf, 2026-10-24, "lanes only where leads can meet"): Coil::terminal_leads_can_meet judges two
// leads by the envelopes of their drawn paths -- polyline, in-section bend arcs, the attach turn's run.
TEST_CASE("terminal_leads_can_meet: envelopes of the drawn paths, bends and turn runs included (ABT #1705)",
          "[constructive-model][coil][abt-1705]") {
    const double od = 1e-3;
    auto radial = [&](double row, int64_t parallel) {
        ConnectionRoute r;
        r.winding = "W";
        r.parallel = parallel;
        r.kind = ConnectionKind::TERMINAL_EXIT;
        r.side = 0;
        r.waypoints = {{10e-3, row}, {20e-3, row}};
        r.plannedBendRadius = 2e-3;
        return r;
    };
    auto at = [](const ConnectionRoute& r) { return std::pair<double, double>{r.waypoints.front().at(1), r.waypoints.front().at(1)}; };
    SECTION("radial leads one coated diameter apart touch, a hair more apart do not") {
        const auto a = radial(0.0, 0);
        const auto touching = radial(od, 1);
        const auto apart = radial(od + 1e-6, 1);
        CHECK(OpenMagnetics::Coil::terminal_leads_can_meet(a, od, at(a), touching, od, at(touching)));
        CHECK_FALSE(OpenMagnetics::Coil::terminal_leads_can_meet(a, od, at(a), apart, od, at(apart)));
    }
    SECTION("a bend cuts its corner: a lead inside the corner meets the arc, not the polyline") {
        // A climbs 10 mm from (10, 0) then runs out radially at 10 mm: a 90 degree corner at (10, 10) with
        // a 4 mm bend, whose arc passes 4 (sqrt(2) - 1) = 1.657 mm inside the corner along its diagonal. B is
        // a point-like lead 2 mm inside the corner on that diagonal: 1.414 mm from both legs (more than one
        // OD) yet 0.343 mm from the arc.
        ConnectionRoute a = radial(0.0, 0);
        a.waypoints = {{10e-3, 0.0}, {10e-3, 10e-3}, {20e-3, 10e-3}};
        a.plannedBendRadius = 4e-3;
        const double inside = 2e-3 / std::sqrt(2.0);
        ConnectionRoute b = radial(0.0, 1);
        b.waypoints = {{10e-3 + inside, 10e-3 - inside}, {10e-3 + inside + 1e-6, 10e-3 - inside}};
        b.plannedBendRadius = 0.5e-3;
        // Without the arc the two stand more than one OD apart.
        ConnectionRoute sharp = a;
        sharp.plannedBendRadius = 0;
        CHECK_FALSE(OpenMagnetics::Coil::terminal_leads_can_meet(sharp, od, at(sharp), b, od, at(b)));
        CHECK(OpenMagnetics::Coil::terminal_leads_can_meet(a, od, at(a), b, od, at(b)));
    }
    SECTION("a steep attach turn climbs over a sibling's lead near the crossing") {
        const auto low = radial(0.0, 0);
        const auto high = radial(1.5e-3, 1);
        CHECK_FALSE(OpenMagnetics::Coil::terminal_leads_can_meet(low, od, at(low), high, od, at(high)));
        // Low's turn climbs to 0.6 mm within the reach: 0.9 mm from high's lead.
        CHECK(OpenMagnetics::Coil::terminal_leads_can_meet(low, od, {0.0, 0.6e-3}, high, od, at(high)));
    }
}

TEST_CASE("terminal_exit_slots: siblings that cannot meet share the crossing, siblings that can stand a pitch apart (ABT #1705)",
          "[constructive-model][coil][abt-1705]") {
    // Alf, 2026-10-24: "lanes only where leads can meet". Section-local leads leaving at heights far
    // apart share x = 0 and leave radially; leads that touch take the side-by-side pitch.
    const double od = 0.855e-3;
    const double radius = 21.4375e-3;
    auto bundle = [&](const std::vector<double>& rows) {
        std::vector<ConnectionRoute> routes;
        for (size_t parallel = 0; parallel < rows.size(); ++parallel) {
            ConnectionRoute r;
            r.winding = "Secondary";
            r.parallel = int64_t(parallel);
            r.kind = ConnectionKind::TERMINAL_EXIT;
            r.side = 0;
            r.waypoints = {{radius, rows[parallel]}, {34.5025e-3, rows[parallel]}};
            r.plannedBendRadius = 0.8275e-3;
            r.fromTurn = "t" + std::to_string(parallel);
            routes.push_back(r);
        }
        return routes;
    };
    const WindingSides secondaryOnly = {{"Secondary", IsolationSide::PRIMARY}};
    SECTION("stacked sections: rows 5 mm apart, gentle advance -> both on the crossing") {
        const std::vector<double> rows = {-20e-3, -15e-3};
        const auto routes = bundle(rows);
        const auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, {od, od}, rows, true, {0.9e-3, 0.9e-3}, round_paths(routes, {od, od}), secondaryOnly);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == 0.0);
    }
    SECTION("touching rows: one pitch apart") {
        const std::vector<double> rows = {-20e-3, -20e-3 + od};
        const auto routes = bundle(rows);
        const auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, {od, od}, rows, true, {0.9e-3, 0.9e-3}, round_paths(routes, {od, od}), secondaryOnly);
        CHECK(std::abs(std::abs(slots[0].value() - slots[1].value()) - od) < 1e-12);
    }
    SECTION("ideal winding keeps the pitch whatever the rows") {
        const std::vector<double> rows = {-20e-3, -15e-3};
        const auto routes = bundle(rows);
        const auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, {od, od}, rows, false, {0.9e-3, 0.9e-3}, {}, secondaryOnly);
        CHECK(std::abs(std::abs(slots[0].value() - slots[1].value()) - od) < 1e-12);
    }
}

TEST_CASE("Real winding: a terminal lead may take the other x half to fit, but a Z face changes only for a different isolation group and only inside the core's exit window (ABT #2274)",
          "[constructive-model][coil][real-winding][abt2274-no-fit-half]") {
    // ABT #2274 (Alf, 2026-11-07, confirmed rule): (1) X: a terminal MAY move to the other x half (-x <-> +x) to
    // avoid other terminals (spec section 22, overflow onto the spare half). (2) Z face: never because of fit
    // alone; a face move only for a group on a DIFFERENT isolation side from the face's other groups (a), an
    // inductor's entrance/exit (b) or the user's pins/sides (c), and only inside the core's exit window. Designs:
    // unmodified MAS examples. 02_flyback_efd25_3c95 and 17_cllc_xfmr_e5528_3c92a: the Primary exit's lane on the
    // primary's half (x <= 0) puts its last turn above the winding window (lane-height check), so it overflows
    // onto the Secondary's half (x >= 0) of the same face. 04_forward_xfmr_e3216_n87: a 3.3V (Tertiary) entrance
    // overflows onto the primary's half of the same face, so no face move. two_switch_forward_transformer_complete:
    // the foil Secondary's axial leads cannot take another lane, do not fit the Primary's face, and that group --
    // a different isolation side -- moves to the other face (section 23), inside the core's exit window. Expected values are derived here from the routes, the windings'
    // isolation sides and the core's columns.
    const std::string design = GENERATE(std::string("02_flyback_efd25_3c95"), std::string("17_cllc_xfmr_e5528_3c92a"),
                                        std::string("04_forward_xfmr_e3216_n87"),
                                        std::string("complete/two_switch_forward_transformer_complete"));
    INFO(design);
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("..")
                           .append("MAS")
                           .append("examples")
                           .append(design + ".json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coilIn);
    auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto layout = coil.get_connection_layout();
    auto group_of = [&](const OpenMagnetics::ConnectionRoute& route) {
        return coil.get_functional_description()[coil.get_winding_index_by_name(route.winding)].get_isolation_side();
    };
    std::set<MAS::IsolationSide> groups;
    std::map<int, std::set<MAS::IsolationSide>> groupsOnFace;
    std::optional<int> primaryFace;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        INFO(label(route));
        REQUIRE(route.pinName.empty());   // the premise: no pin decides a half, so S0/S1 alone do
        groups.insert(group_of(route));
        groupsOnFace[route.side].insert(group_of(route));
        if (group_of(route) == MAS::IsolationSide::PRIMARY) {
            REQUIRE((!primaryFace || primaryFace.value() == route.side));
            primaryFace = route.side;
        }
    }
    REQUIRE(groups.size() > 1);
    REQUIRE(primaryFace.has_value());
    // Rule 2a: one isolation side never parts faces for fit, and the Primary's group (S0's first) keeps S0's face.
    for (const auto& side : groups) {
        size_t faces = 0;
        for (const auto& [face, onFace] : groupsOnFace) {
            faces += onFace.count(side);
        }
        INFO("isolation side " << magic_enum::enum_name(side));
        CHECK(faces == 1);
    }
    CHECK(primaryFace.value() == 0);
    // Rule 1: the half is S0's (the group's) on a face shared by groups, S1's (entrance x <= 0, exit x >= 0) on a
    // face one group has to itself (spec section 24). A lead that overflowed (spec section 22) stands on the OTHER
    // half of the SAME face (an x move, never a face move); every other lead stands on its own half.
    size_t judged = 0;
    size_t overflowed = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        INFO(label(route) << " face " << route.side);
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        const int expected = groupsOnFace.at(route.side).size() == 1 ? (entrance ? -1 : 1)
                                                                      : (group_of(route) == MAS::IsolationSide::PRIMARY ? -1 : 1);
        REQUIRE(route.laneHalf.has_value());
        CHECK(route.laneHalf.value() == expected);
        REQUIRE(route.laneOverflow.has_value());
        if (route.exitSlot && !route.axialLead) {
            INFO("lane x = " << route.exitSlot.value() * 1e3 << " mm, overflow " << route.laneOverflow.value());
            if (route.laneOverflow.value()) {
                CHECK(expected * route.exitSlot.value() < -1e-9);
                ++overflowed;
            }
            else {
                CHECK(expected * route.exitSlot.value() >= -1e-9);
            }
            ++judged;
        }
        else {
            CHECK_FALSE(route.laneOverflow.value());
        }
    }
    CHECK(judged >= 2);
    if (design.starts_with("02_") || design.starts_with("17_")) {
        // The Primary exit does not fit the primary's half (its Z dragback holds the crossing, and every lane
        // further out lifts its last turn above the window): it overflows onto the Secondary's half of the SAME
        // face and fits there. Nothing moved face: every group stays on S0's face.
        CHECK(groupsOnFace.size() == 1);
        CHECK(overflowed >= 1);
        bool found = false;
        for (const auto& route : layout.routes) {
            if (is_terminal(route) && route.winding == "Primary" && route.kind == ConnectionKind::TERMINAL_EXIT) {
                found = true;
                INFO(label(route));
                REQUIRE(route.exitSlot.has_value());
                CHECK(route.laneOverflow == std::optional<bool>(true));
                CHECK(route.exitSlot.value() > 1e-9);
                CHECK(route.side == primaryFace.value());
                CHECK_FALSE(std::any_of(route.findings.begin(), route.findings.end(),
                                        [](const auto& finding) { return finding.kind == LeadFindingKind::LANES_DO_NOT_FIT; }));
            }
        }
        CHECK(found);
    }
    else if (design.starts_with("04_")) {
        // The 3.3V (Tertiary) leads do not all fit the non-primary half of the shared face: the one that does not
        // overflows onto the primary's half of the SAME face (rule 1, an x move) and fits there, so no group
        // changes face for it (rule 2: the x move is tried first; a face move is not taken for fit alone).
        CHECK(groupsOnFace.size() == 1);
        CHECK(overflowed >= 1);
        for (const auto& route : layout.routes) {
            if (!is_terminal(route)) {
                continue;
            }
            INFO(label(route));
            if (route.laneOverflow == std::optional<bool>(true)) {
                CHECK(group_of(route) == MAS::IsolationSide::TERTIARY);
            }
            for (const auto& finding : route.findings) {
                CHECK(finding.kind != LeadFindingKind::LANES_DO_NOT_FIT);
            }
        }
    }
    else {
        // two_switch: the foil Secondary's axial leads (which no x lane can move) do not fit the Primary's face.
        // Rule 2: the group that did not fit the Primary's face is a DIFFERENT isolation side from it (a); it moved
        // to the other face, and it fits there: no misfit is left, and every moved lead's radial run, from the
        // bobbin column's face outwards, stays clear of the core's outer legs (E cores: the lateral columns, which
        // span the core's depth) -- the core's exit window.
        const auto moved = MAS::IsolationSide::SECONDARY;
        const auto bobbin = coil.resolve_bobbin();
        REQUIRE(bobbin.get_processed_description());
        const double columnHalfDepth = bobbin.get_processed_description()->get_column_depth();
        const auto columns = core.get_columns();
        size_t movedLeads = 0;
        for (const auto& route : layout.routes) {
            if (!is_terminal(route)) {
                continue;
            }
            INFO(label(route) << " face " << route.side);
            for (const auto& finding : route.findings) {
                CHECK(finding.kind != LeadFindingKind::LANES_DO_NOT_FIT);
            }
            if (group_of(route) != moved) {
                continue;
            }
            CHECK(route.side != primaryFace.value());
            CHECK(group_of(route) != MAS::IsolationSide::PRIMARY);
            CHECK(groupsOnFace.at(route.side).count(MAS::IsolationSide::PRIMARY) == 0);
            ++movedLeads;
            const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
            const double lane = route.axialLead ? (entrance ? -1.0 : 1.0) * route.axialLead->alongFace : route.exitSlot.value();
            double radius = 0;
            if (route.axialLead) {
                radius = route.axialLead->outerDiameter / 2;
            }
            else {
                const auto turns = coil.get_turns_description().value();
                const auto turn = std::find_if(turns.begin(), turns.end(), [&](const auto& candidate) {
                    return candidate.get_name() == (entrance ? route.toTurn : route.fromTurn);
                });
                REQUIRE(turn != turns.end());
                radius = coil.terminal_turn_path(route, *turn, layout.routes).stubDiameter / 2;
                REQUIRE(radius > 0);
            }
            for (size_t c = 1; c < columns.size(); ++c) {
                REQUIRE(columns[c].get_shape() == MAS::ColumnShape::RECTANGULAR);   // E-type legs: boxes
                const auto at = columns[c].get_coordinates();
                const double cz = at.size() > 2 ? at[2] : 0.0;
                const bool overlapsX = std::abs(lane - at[0]) < columns[c].get_width() / 2 + radius;
                const bool reachesRun = std::abs(cz) + columns[c].get_depth() / 2 > columnHalfDepth + 1e-12;
                CHECK_FALSE((overlapsX && reachesRun));
            }
        }
        CHECK(movedLeads >= 2);
    }
    settings.reset();
}

// AGENT11-SCAN (temporary, remove before commit): faces/halves/findings of every MAS example's terminal leads.
TEST_CASE("agent11 scan of terminal faces", "[.][agent11-scan]") {
    namespace fs = std::filesystem;
    const auto root = fs::path{std::source_location::current().file_name()}.parent_path().append("..").append("MAS").append("examples");
    std::vector<fs::path> files;
    for (const auto& dir : {root, fs::path(root).append("complete")}) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() == ".json") {
                files.push_back(entry.path());
            }
        }
    }
    std::sort(files.begin(), files.end());
    const char* only = std::getenv("AGENT11_SCAN_ONLY");
    for (const auto& path : files) {
        if (only && path.stem().string().find(only) == std::string::npos) {
            continue;
        }
        std::cout << "=== " << path.stem().string() << std::endl;
        try {
            std::ifstream file(path);
            const auto masJson = json::parse(file);
            const auto& magneticJson = masJson.at("magnetic");
            auto& settings = Settings::GetInstance();
            settings.reset();
            settings.set_coil_use_real_winding_geometry(true);
            settings.set_coil_lead_bend_radius_factor(1.05);
            settings.set_coil_lead_minimum_bend_radius(std::nullopt);
            OpenMagnetics::Core core(magneticJson.at("core"));
            OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coilIn);
            auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
            std::cout << "  family " << magic_enum::enum_name(core.get_shape_family()) << " blocking " << coil.is_real_winding_blocking_applied() << std::endl;
            for (const auto& route : coil.get_connection_layout().routes) {
                if (!is_terminal(route)) {
                    continue;
                }
                const auto side = coil.get_functional_description()[coil.get_winding_index_by_name(route.winding)].get_isolation_side();
                std::cout << "  " << label(route) << " iso " << magic_enum::enum_name(side) << " face " << route.side << " slot "
                          << (route.exitSlot ? std::to_string(route.exitSlot.value() * 1e3) : std::string("none")) << " half "
                          << (route.laneHalf ? std::to_string(route.laneHalf.value()) : std::string("none")) << " overflow "
                          << (route.laneOverflow ? std::to_string(route.laneOverflow.value()) : std::string("none"))
                          << (route.axialLead ? " axial" : "") << std::endl;
                for (const auto& finding : route.findings) {
                    std::cout << "    " << magic_enum::enum_name(finding.kind) << ": " << finding.witness << std::endl;
                }
            }
        }
        catch (const std::exception& error) {
            std::cout << "  THROWS: " << error.what() << std::endl;
        }
    }
    Settings::GetInstance().reset();
}

// ABT #2274 (Alf, 2026-11-07, "slot-step ramp"): a terminal lead whose stub is not drawable (its route
// starts at the row) still leaves its turn at its settled slot, where the turn stands off the row by its
// uniform climb. MKF plans that step as a cosine S on the turn (ConnectionRoute::slotRampLength /
// slotRampHeight) without moving a lane. Design: the UNMODIFIED MAS example 02_flyback_efd25_3c95; the
// steps are the ones MVB++ drew as a slanted elbow leg at the agent-12 state (agent 13's measurements):
// Primary exit 117 um, Secondary p2 exit 35 um, Secondary p1 entrance and exit 17.7 um, Primary entrance
// 4.66 um. Its dragbacks: each bend at the lead bend radius of the wire, the step-out never past one OD.
TEST_CASE("Real winding: a non-drawable terminal stub climbs its slot step as a ramp, and a dragback's bump is one OD (ABT #2274)",
          "[constructive-model][coil][real-winding][abt2274-slot-ramp]") {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("..")
                           .append("MAS")
                           .append("examples")
                           .append("02_flyback_efd25_3c95.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coilIn(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coilIn);
    auto coil = OpenMagnetics::magnetic_autocomplete(magnetic, json{}).get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto layout = coil.get_connection_layout();
    const auto windings = coil.get_functional_description();
    const auto turnsAll = coil.get_turns_description().value();
    auto find_route = [&](size_t windingIndex, int64_t parallel, OpenMagnetics::ConnectionKind kind) -> const OpenMagnetics::ConnectionRoute* {
        for (const auto& route : layout.routes) {
            if (route.winding == windings.at(windingIndex).get_name() && route.parallel == parallel && route.kind == kind) {
                return &route;
            }
        }
        return nullptr;
    };
    for (const auto& route : layout.routes) {
        if (route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE || route.kind == OpenMagnetics::ConnectionKind::TERMINAL_EXIT) {
            std::cout << "[abt2274-slot-ramp] " << route.winding << " p" << route.parallel
                      << (route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE ? " IN" : " OUT") << " waypoints "
                      << route.waypoints.size() << " slot " << route.exitSlot.value_or(std::nan("")) * 1e3 << " mm height "
                      << route.slotRampHeight.value_or(std::nan("")) * 1e6 << " um length "
                      << route.slotRampLength.value_or(std::nan("")) * 1e3 << " mm lift " << route.liftLength.has_value()
                      << " axial " << route.axialLead.has_value() << " tangential " << route.laneTangential.value_or(false)
                      << " side " << route.side << " bend " << route.plannedBendRadius * 1e3 << " row "
                      << route.waypoints.front().at(1) * 1e3 << " turn " << (route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE ? route.toTurn : route.fromTurn)
                      << " findings " << route.findings.size() << std::endl;
            for (const auto& finding : route.findings) {
                std::cout << "    finding " << static_cast<int>(finding.kind) << ": " << finding.witness.substr(0, 400) << std::endl;
            }
        }
    }
    struct Expected {
        size_t winding;
        int64_t parallel;
        OpenMagnetics::ConnectionKind kind;
        double step;   // [m], magnitude
    };
    const std::vector<Expected> expected{
        // The leads agent 13 found stepping at their slots (MVB++ drew 117 / 35 / 17.7 / 17.7 / 4.66 um of
        // slant in its 1-OD elbow leg, a part of the step); the steps themselves, row - turn height at
        // the slot, are pinned here and re-derived independently below from the turns description.
        {0, 0, OpenMagnetics::ConnectionKind::TERMINAL_EXIT, 144.552e-6},
        {1, 2, OpenMagnetics::ConnectionKind::TERMINAL_EXIT, 47.5283e-6},
        {1, 1, OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE, 23.7641e-6},
        {1, 1, OpenMagnetics::ConnectionKind::TERMINAL_EXIT, 23.7641e-6},
        {0, 0, OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE, 5.29802e-6},
    };
    for (const auto& lead : expected) {
        const auto* route = find_route(lead.winding, lead.parallel, lead.kind);
        REQUIRE(route != nullptr);
        INFO(route->winding << " p" << route->parallel << (lead.kind == OpenMagnetics::ConnectionKind::TERMINAL_EXIT ? " exit" : " entrance"));
        // The premise: a stub that is not drawable, so its route starts at the row.
        REQUIRE(route->waypoints.size() == 2);
        REQUIRE(route->exitSlot);
        REQUIRE(!route->rampLength);
        REQUIRE(route->slotRampHeight);
        REQUIRE(route->slotRampLength);
        CHECK(std::abs(route->slotRampHeight.value()) == Catch::Approx(lead.step).epsilon(1e-4));
        // Independently from the turns description: the attach turn's station and its advance to the
        // conductor's neighbouring turn, climbed uniformly over the turn's length; on the flat face the
        // lane x is the distance along the turn from the crossing (upstream for an entrance).
        {
            const bool entrance = lead.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE;
            const std::string turnName = entrance ? route->toTurn : route->fromTurn;
            std::optional<size_t> at;
            for (size_t t = 0; t < turnsAll.size(); ++t) {
                if (turnsAll[t].get_name() == turnName) {
                    at = t;
                }
            }
            REQUIRE(at);
            std::optional<double> neighbourY;
            if (entrance) {
                for (size_t t = at.value() + 1; t < turnsAll.size() && !neighbourY; ++t) {
                    if (turnsAll[t].get_winding() == turnsAll[at.value()].get_winding() && turnsAll[t].get_parallel() == turnsAll[at.value()].get_parallel()) {
                        neighbourY = turnsAll[t].get_coordinates().at(1);
                    }
                }
            }
            else {
                for (size_t t = at.value(); t-- > 0 && !neighbourY;) {
                    if (turnsAll[t].get_winding() == turnsAll[at.value()].get_winding() && turnsAll[t].get_parallel() == turnsAll[at.value()].get_parallel()) {
                        neighbourY = turnsAll[t].get_coordinates().at(1);
                    }
                }
            }
            REQUIRE(neighbourY);
            const double station = turnsAll[at.value()].get_coordinates().at(1);
            const double advance = entrance ? neighbourY.value() - station : station - neighbourY.value();
            // A conductor's first turn may state no length of its own; its neighbour's is the same turn length.
            double turnLength = turnsAll[at.value()].get_length();
            for (size_t t = 0; t < turnsAll.size() && !(turnLength > 0); ++t) {
                if (turnsAll[t].get_winding() == turnsAll[at.value()].get_winding() && turnsAll[t].get_parallel() == turnsAll[at.value()].get_parallel() &&
                    turnsAll[t].get_coordinates().at(1) == neighbourY.value()) {
                    turnLength = turnsAll[t].get_length();
                }
            }
            REQUIRE(turnLength > 0);
            const double along = route->exitSlot.value();
            const double independent = route->waypoints.front().at(1) - (station - advance * along / turnLength);
            std::cout << "[abt2274-slot-ramp] independent step " << independent * 1e6 << " um vs MKF " << route->slotRampHeight.value() * 1e6
                      << " um (advance " << advance * 1e6 << " um over " << turnLength * 1e3 << " mm)" << std::endl;
            // The turn length here is the wound turn's own; MKF's face path adds the ride over the other face
            // (TerminalTurnPath::rideLength), ~2 % on 02, so the two agree to that.
            CHECK(route->slotRampHeight.value() == Catch::Approx(independent).epsilon(0.03));
        }
        // The S is the ramp rule's own: S = pi sqrt(h Rg / 2), Rg the plannedBendRadius on a flat face
        // (lead_ramp_length), so its curvature peaks at exactly 1 / Rg.
        const double h = std::abs(route->slotRampHeight.value());
        const double S = route->slotRampLength.value();
        CHECK(h * std::numbers::pi * std::numbers::pi / (2 * S * S) == Catch::Approx(1.0 / route->plannedBendRadius).epsilon(1e-9));
        // ABT #2274 (agent 17): MKF publishes where the S ends (the consumer draws onto it): bend + S from the lane,
        // on the side the turn covers.
        REQUIRE(route->slotRampEndX);
        REQUIRE(route->slotRampEndHeight);
        CHECK(std::abs(route->slotRampEndX.value() - route->exitSlot.value()) == Catch::Approx(S + route->plannedBendRadius).epsilon(1e-12));
    }
    // Dragbacks: bent at the wire's lead bend radius, the step-out stopping exactly at the return's radius.
    size_t dragbacks = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != OpenMagnetics::ConnectionKind::Z_DRAGBACK) {
            continue;
        }
        ++dragbacks;
        INFO(route.winding << " p" << route.parallel << " dragback " << route.fromTurn << " -> " << route.toTurn);
        REQUIRE(route.dragbackBendRadius);
        REQUIRE(route.dragbackBendAngle);
        const size_t windingIndex = coil.get_winding_index_by_name(route.winding);
        const auto wire = coil.resolve_wire(windingIndex);
        const double od = std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height());
        const double R = route.dragbackBendRadius.value();
        const double theta = route.dragbackBendAngle.value();
        const double step = std::abs(route.waypoints.back().at(0) - route.waypoints.front().at(0));
        // 02's dragbacks step one OD out (adjacent layers), and bend at the radius a lead of the same
        // wire bends at (its unsleeved terminal lead's plannedBendRadius), never below the declared 1.05 x.
        CHECK(step == Catch::Approx(od).epsilon(1e-6));
        CHECK(R >= 1.05 * od / 2);
        const auto* lead = find_route(windingIndex, route.parallel, OpenMagnetics::ConnectionKind::TERMINAL_EXIT);
        REQUIRE(lead != nullptr);
        REQUIRE(!lead->sleeveOuterDiameter);
        CHECK(R == Catch::Approx(lead->plannedBendRadius).epsilon(1e-12));
        // The bump: two bends of theta reach exactly the step, never past it.
        const double reach = theta < std::numbers::pi / 2 ? 2 * R * (1 - std::cos(theta)) : step;
        CHECK(reach == Catch::Approx(step).epsilon(1e-9));
        CHECK(theta <= std::numbers::pi / 2);
    }
    CHECK(dragbacks > 0);
    Settings::GetInstance().reset();
}

// ABT #2274 (Alf, 2026-11-07: "replace the 65-point sampling with a RIGOROUS bound"): the slot-ramp fit
// check's distance, slot_ramp_closest_approach, against a dense brute force, and on a thin-wire pair whose
// closest approach falls between the 65 stations the earlier check sampled (that sampling reported it clear).
TEST_CASE("Real winding: a slot ramp's closest approach to another lead is bounded rigorously (ABT #2274)",
          "[constructive-model][coil][real-winding][abt2274-slot-ramp]") {
    using OpenMagnetics::SlotRampCentreline;
    const SlotRampCentreline ramp{0.0, 1.0, 5e-3, 1e-3, 144.552e-6, 0.544e-3, 0.415e-3};
    auto point_on = [](const SlotRampCentreline& c, double s) { return OpenMagnetics::slot_ramp_point(c, s); };
    auto point_segment = [](double r, double y, const std::vector<double>& a, const std::vector<double>& b) {
        const double ur = b[0] - a[0], uy = b[1] - a[1];
        const double t = std::clamp(((r - a[0]) * ur + (y - a[1]) * uy) / (ur * ur + uy * uy), 0.0, 1.0);
        return std::hypot(r - a[0] - t * ur, y - a[1] - t * uy);
    };
    auto brute_route = [&](double otherLaneX, const std::vector<std::vector<double>>& route, int stations) {
        double best = std::numeric_limits<double>::max();
        const double end = ramp.rampLength + ramp.bendRadius;
        for (int k = 0; k <= stations; ++k) {
            const auto p = point_on(ramp, end * k / stations);
            double section = std::numeric_limits<double>::max();
            for (size_t w = 0; w + 1 < route.size(); ++w) {
                section = std::min(section, point_segment(p[1], p[2], route[w], route[w + 1]));
            }
            best = std::min(best, std::hypot(p[0] - otherLaneX, section));
        }
        return best;
    };
    SECTION("against a route: matches a dense brute force, fires just under it, not just over it") {
        const std::vector<std::vector<double>> route{{5.3e-3, 0.95e-3}, {5.3e-3, 2e-3}, {7e-3, 2e-3}};
        const double otherLaneX = 0.3e-3;
        const double dense = brute_route(otherLaneX, route, 400000);
        const auto found = OpenMagnetics::slot_ramp_closest_approach(ramp, otherLaneX, route, std::nullopt, dense + 1e-6);
        CHECK(found.distance < dense + 1e-6);
        CHECK(found.distance >= dense - 3e-9);   // never under the true minimum (the brute force overestimates by <= 1.5 nm)
        CHECK(found.otherAlong < 0);
        const auto clear = OpenMagnetics::slot_ramp_closest_approach(ramp, otherLaneX, route, std::nullopt, dense - 1e-6);
        CHECK_FALSE(clear.distance < dense - 1e-6);
    }
    SECTION("a thin-wire touch between two of the 65 stations is found") {
        const double sStar = ramp.bendRadius + ramp.rampLength * 0.5 / 64;   // halfway between the S's stations 0 and 1
        const auto p = point_on(ramp, sStar);
        const double gap = 18e-6, clearance = 18.4e-6;   // two 0.0184 mm coated wires, 0.4 um overlap at s*
        const std::vector<std::vector<double>> route{{p[1] + gap, p[2] - 1e-9}, {p[1] + gap, p[2] + 1e-9}};
        double sampled = std::numeric_limits<double>::max();
        for (int k = 0; k <= 64; ++k) {
            const auto q = point_on(ramp, ramp.bendRadius + ramp.rampLength * k / 64.0);
            sampled = std::min(sampled, std::hypot(q[0] - p[0], point_segment(q[1], q[2], route[0], route[1])));
        }
        INFO("65-station sampling: " << sampled * 1e6 << " um");
        CHECK_FALSE(sampled < clearance);
        const auto found = OpenMagnetics::slot_ramp_closest_approach(ramp, p[0], route, std::nullopt, clearance);
        CHECK(found.distance < clearance);
        CHECK(found.distance >= gap - 1e-12);
        // the exhibited point is within the overlap's reach of s*: |s - s*| < sqrt(clearance^2 - gap^2)
        CHECK(std::abs(found.along - sStar) < std::sqrt(clearance * clearance - gap * gap));
    }
    SECTION("against the other lead's own slot ramp: matches a dense 2-D brute force") {
        const SlotRampCentreline other{1.2e-3, -1.0, 5.25e-3, 1.1e-3, -60e-6, 0.4e-3, 0.3e-3};
        const std::vector<std::vector<double>> route{{5.25e-3, 1.1e-3}, {5.25e-3, 3e-3}};
        const int n = 1500;
        double dense = brute_route(other.laneX, route, 20000);
        for (int i = 0; i <= n; ++i) {
            const auto a = point_on(ramp, (ramp.rampLength + ramp.bendRadius) * i / n);
            for (int j = 0; j <= n; ++j) {
                const auto b = point_on(other, (other.rampLength + other.bendRadius) * j / n);
                dense = std::min(dense, std::hypot(a[0] - b[0], std::hypot(a[1] - b[1], a[2] - b[2])));
            }
        }
        const auto found = OpenMagnetics::slot_ramp_closest_approach(ramp, other.laneX, route, other, dense + 1e-6);
        CHECK(found.distance < dense + 1e-6);
        CHECK(found.distance > dense - 2e-6);   // the grid brute force overestimates by at most its station spacing
        CHECK(found.otherAlong >= 0);
    }
}
