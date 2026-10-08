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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <source_location>

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

    size_t crossBundlePairs = 0;
    size_t primaryEntranceVsExit = 0;
    double minimumCrossBundleSurplus = std::numeric_limits<double>::max();
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
            // Two bundles whose routes come within reach of each other's corners (the coated radii plus
            // both bends) stand the coated radii plus the larger bend apart along the face: disjoint,
            // with a clearance gap, never an exact touch.
            if (section < radii + a.plannedBendRadius + b.plannedBendRadius) {
                ++crossBundlePairs;
                if (a.winding == "Primary" && b.winding == "Primary") {
                    ++primaryEntranceVsExit;
                }
                const double required = radii + std::max(a.plannedBendRadius, b.plannedBendRadius);
                CHECK(gap >= required - 1e-9);
                CHECK(gap > radii);
                minimumCrossBundleSurplus = std::min(minimumCrossBundleSurplus, gap - radii);
            }
        }
    }
    // The pairs this test is about exist: the Primary's inner-layer entrances against its outer-layer exits.
    REQUIRE(primaryEntranceVsExit == 16);
    REQUIRE(crossBundlePairs >= primaryEntranceVsExit);
    INFO("smallest cross-bundle clearance beyond the coated radii: " << minimumCrossBundleSurplus * 1e3 << " mm");
    CHECK(minimumCrossBundleSurplus > 0);

    // Each lead leaves its helix with room for its bends: a straight stub from its turn to its row
    // holds two bends of plannedBendRadius (one in the turn surface, one onto the radial run), so it is
    // at least 2R long, or it is published as a ramp on the turn surface (ABT #1336) whose span the
    // lanes above already keep clear.
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
        const double height = std::abs(turnEnd[1] - next[1]);
        INFO(label(lead) << ": stub " << height * 1e3 << " mm, bend radius " << lead.plannedBendRadius * 1e3 << " mm");
        CHECK((height >= 2 * lead.plannedBendRadius - 1e-12 || lead.rampLength.has_value()));
    }
    REQUIRE(stubs > 0);
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
    // MKF's own lanes pass.
    CHECK_NOTHROW(OpenMagnetics::Coil::check_terminal_lanes(routes, diameters, slots, true));
    // The lanes the old rule handed out put the Primary's parallel 3 exit on its parallel 0 entrance's lane:
    // their rows touch in the window section (0.855 mm), so the old test saw no clash.
    slots[exit.value()] = slots[entrance.value()];
    CHECK_THROWS_WITH(OpenMagnetics::Coil::check_terminal_lanes(routes, diameters, slots, true),
                      Catch::Matchers::ContainsSubstring("Terminal lanes collide") &&
                          Catch::Matchers::ContainsSubstring("Primary' parallel 3 exit") &&
                          Catch::Matchers::ContainsSubstring("different bundles"));
    // Two leads on one point collide even without the bend rule.
    auto coincident = routes;
    coincident[exit.value()].waypoints = coincident[entrance.value()].waypoints;
    CHECK_THROWS_WITH(OpenMagnetics::Coil::check_terminal_lanes(coincident, diameters, slots, false),
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

    SECTION("ideal winding: the rows touch, so the old rule lets both bundles start at the plane") {
        auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, false);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == od);
        CHECK(slots[3].value() == 0.0);   // exit order by span descending: parallel 1 first
        CHECK(slots[2].value() == od);
    }
    SECTION("real winding: the exit bundle starts the first lane at least d + R past the entrance bundle") {
        auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, true);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == od);
        // The second entrance stands at od, so the first exit needs od + (od + bend) = 2.159 mm; on the pitch grid
        // that is lane 3 = 2.565 mm.
        const double firstExitLane = std::ceil((od + od + bend) / od - 1e-9) * od;
        CHECK(firstExitLane == 3 * od);
        CHECK(slots[3].value() == firstExitLane);
        CHECK(slots[2].value() == firstExitLane + od);
    }
    SECTION("real winding without a planned bend throws") {
        routes[2].plannedBendRadius = 0;
        CHECK_THROWS_WITH(OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, true),
                          Catch::Matchers::ContainsSubstring("planned bend radius"));
    }
}
